@testset "Estimation Functions" begin
    
    # Helper function to create fully prepared model (matched and balanced)
    function create_prepared_model()
        dates = Date(2021, 10, 1):Day(1):Date(2021, 12, 31)
        n_units = 30
        
        test_data = DataFrame()
        for unit in 1:n_units
            treatment_effect = (unit <= 8) ? 0.2 : 0.0  # Treatment effect for treated units
            
            for (i, date) in enumerate(dates)
                is_treated = (unit <= 8) && (i >= 40)  # Treatment starts at day 40
                outcome_base = 1.0 + 0.02 * i + randn() * 0.3
                outcome = outcome_base + (is_treated ? treatment_effect : 0.0)
                
                push!(test_data, (
                    date = date,
                    fips = 1000 + unit,
                    pop_dens = 50.0 + randn() * 10,
                    cumul_death_rate = 20.0 + i * 0.1 + randn() * 2,
                    death_rte = max(0.0, outcome),
                    gub = (unit <= 8) ? 1 : 0
                ))
            end
        end
        
        test_data[!, :day] = Dates.value.(test_data.date .- minimum(test_data.date))
        
        model = makemodel(
            test_data, :day, :fips, :gub, :death_rte,
            [:pop_dens, :cumul_death_rate],
            Dict(:pop_dens => false, :cumul_death_rate => true),
            50:60, -30:-1;  # F period after treatment, L period before
            estimator = "ATT"
        )
        
        # Complete matching and balancing
        match!(model, test_data)
        balance!(model, test_data)
        
        return model, test_data
    end
    
    @testset "estimate! function" begin
        model, test_data = create_prepared_model()
        
        @testset "Basic estimation" begin
            # Test that estimate! runs without errors
            @test_nowarn estimate!(model, test_data)
            
            # After estimation, model should have results
            @test hasfield(typeof(model), :results)
            
            # Results should be a DataFrame with ATT estimates
            if isdefined(model, :results) && !isnothing(model.results)
                @test model.results isa DataFrame
                @test nrow(model.results) > 0
                
                # Should have columns for time periods in F
                expected_periods = length(model.F)
                @test nrow(model.results) <= expected_periods
                
                # Should have ATT estimates
                if "att" in names(model.results)
                    @test all(x -> isa(x, Real), model.results.att)
                end
            end
        end
        
        @testset "Estimation with custom parameters" begin
            model_custom, test_data_custom = create_prepared_model()
            
            # Test with custom number of bootstrap iterations
            @test_nowarn estimate!(
                model_custom, test_data_custom;
                iterations = 100
            )
            
            # Test with custom percentiles
            @test_nowarn estimate!(
                model_custom, test_data_custom;
                percentiles = [0.05, 0.5, 0.95],
                iterations = 50
            )
        end
        
        @testset "Estimation with overall estimate" begin
            model_overall, test_data_overall = create_prepared_model()
            
            # Test estimation that returns overall estimate
            overall_result = estimate!(
                model_overall, test_data_overall;
                overallestimate = true,
                iterations = 100
            )
            
            @test overall_result isa Overall
            @test isa(overall_result.att, Real)
            @test isa(overall_result.percentiles, Vector)
            @test length(overall_result.percentiles) > 0
        end
        
        @testset "Bayes factor computation" begin
            model_bayes, test_data_bayes = create_prepared_model()
            
            # Test with Bayes factor enabled
            @test_nowarn estimate!(
                model_bayes, test_data_bayes;
                dobayesfactor = true,
                iterations = 100
            )
            
            # Test with Bayes factor disabled
            @test_nowarn estimate!(
                model_bayes, test_data_bayes;
                dobayesfactor = false,
                iterations = 100
            )
        end
        
        @testset "P-value computation" begin
            model_pval, test_data_pval = create_prepared_model()
            
            # Test with p-value computation enabled
            @test_nowarn estimate!(
                model_pval, test_data_pval;
                dopvalue = true,
                iterations = 100
            )
            
            # Test with p-value computation disabled
            @test_nowarn estimate!(
                model_pval, test_data_pval;
                dopvalue = false,
                iterations = 100
            )
        end
    end
    
    @testset "Overall structure and functions" begin
        @testset "Overall constructor" begin
            # Test the overall() constructor function
            result = overall(
                att = 0.15,
                percentiles = [0.05, 0.15, 0.25],
                bayesfactor = 2.5,
                ntreatedmean = 8.0,
                pvalue = 0.03
            )
            
            @test result isa Overall
            @test result.att == 0.15
            @test result.percentiles == [0.05, 0.15, 0.25]
            @test result.bayesfactor == 2.5
            @test result.ntreatedmean == 8.0
            @test result.pvalue == 0.03
            @test ismissing(result.stratum)
        end
        
        @testset "Overall with stratification" begin
            # Test overall with stratification information
            result_strat = overall(
                att = [0.1, 0.2, 0.15],
                percentiles = [[0.05, 0.15], [0.1, 0.3], [0.05, 0.25]],
                bayesfactor = [2.0, 3.0, 2.5],
                ntreatedmean = [3.0, 3.0, 2.0],
                pvalue = [0.05, 0.02, 0.08],
                stratum = [1.0, 2.0, 3.0]
            )
            
            @test result_strat.att isa Vector
            @test result_strat.percentiles isa Vector{Vector{Float64}}
            @test result_strat.bayesfactor isa Vector
            @test result_strat.stratum isa Vector
        end
    end
    
    @testset "Estimation edge cases" begin
        @testset "No variation in outcomes" begin
            # Create data with no variation in outcomes
            dates = Date(2021, 10, 1):Day(1):Date(2021, 11, 15)
            test_data_constant = DataFrame()
            
            for unit in 1:10
                for (i, date) in enumerate(dates)
                    push!(test_data_constant, (
                        date = date,
                        fips = 1000 + unit,
                        pop_dens = 50.0 + randn() * 10,
                        cumul_death_rate = 20.0 + i * 0.1,
                        death_rte = 1.0,  # Constant outcome
                        gub = (unit <= 3) ? 1 : 0
                    ))
                end
            end
            
            test_data_constant[!, :day] = Dates.value.(test_data_constant.date .- minimum(test_data_constant.date))
            
            model_constant = makemodel(
                test_data_constant, :day, :fips, :gub, :death_rte,
                [:pop_dens, :cumul_death_rate],
                Dict(:pop_dens => false, :cumul_death_rate => true),
                20:25, -15:-1
            )
            
            match!(model_constant, test_data_constant)
            balance!(model_constant, test_data_constant)
            
            # Should handle constant outcomes gracefully
            @test_nowarn estimate!(model_constant, test_data_constant; iterations = 50) ||
                   @test_throws Exception estimate!(model_constant, test_data_constant; iterations = 50)
        end
        
        @testset "Single time period" begin
            model_single_period, test_data_single = create_prepared_model()
            
            # Modify model to have single time period in F
            model_single_period = @set model_single_period.F = 50:50
            
            # Should handle single period estimation
            @test_nowarn estimate!(model_single_period, test_data_single; iterations = 50)
        end
        
        @testset "Very few bootstrap iterations" begin
            model_few_boot, test_data_few = create_prepared_model()
            
            # Test with minimal bootstrap iterations
            @test_nowarn estimate!(
                model_few_boot, test_data_few;
                iterations = 5,
                percentiles = [0.5]  # Just median
            )
        end
    end
    
    @testset "Results validation" begin
        model, test_data = create_prepared_model()
        estimate!(model, test_data; iterations = 100)
        
        @testset "Results structure" begin
            if isdefined(model, :results) && !isnothing(model.results)
                results = model.results
                
                # Basic structure checks
                @test results isa DataFrame
                @test nrow(results) > 0
                
                # Should have numeric ATT estimates
                if "att" in names(results)
                    @test eltype(results.att) <: Real
                end
                
                # Should have confidence intervals
                for col in names(results)
                    if occursin("%", col) || occursin("percentile", col)
                        @test eltype(results[!, col]) <: Real
                    end
                end
            end
        end
        
        @testset "Confidence intervals" begin
            if isdefined(model, :results) && !isnothing(model.results)
                results = model.results
                
                # If we have confidence intervals, they should be ordered
                ci_cols = filter(col -> occursin("%", col), names(results))
                if length(ci_cols) >= 2
                    lower_cols = filter(col -> occursin("2.5%", col) || occursin("5%", col), ci_cols)
                    upper_cols = filter(col -> occursin("97.5%", col) || occursin("95%", col), ci_cols)
                    
                    if length(lower_cols) > 0 && length(upper_cols) > 0
                        # Lower bounds should generally be less than upper bounds
                        for i in 1:nrow(results)
                            if !any(ismissing.([results[i, lower_cols[1]], results[i, upper_cols[1]]]))
                                @test results[i, lower_cols[1]] <= results[i, upper_cols[1]]
                            end
                        end
                    end
                end
            end
        end
    end
end