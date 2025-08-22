@testset "Model Construction" begin
    
    @testset "makemodel function" begin
        # Create test data similar to vignette
        dates = Date(2021, 10, 1):Day(1):Date(2021, 12, 31)
        n_units = 50
        n_days = length(dates)
        
        # Create synthetic test data
        test_data = DataFrame()
        for unit in 1:n_units
            for (i, date) in enumerate(dates)
                push!(test_data, (
                    date = date,
                    fips = 1000 + unit,
                    pop_dens = 50.0 + randn() * 10,
                    cumul_death_rate = 20.0 + i * 0.1 + randn() * 2,
                    death_rte = abs(randn() * 0.5),
                    gub = (unit <= 10) ? 1 : 0  # First 10 units are treated
                ))
            end
        end
        
        # Add integer time variable
        test_data[!, :day] = Dates.value.(test_data.date .- minimum(test_data.date))
        
        # Test basic model construction
        @testset "Basic model creation" begin
            matching_covariates = [:pop_dens, :cumul_death_rate]
            timevary = Dict(:pop_dens => false, :cumul_death_rate => true)
            
            model = makemodel(
                test_data, :day, :fips, :gub, :death_rte,
                matching_covariates, timevary,
                10:20, -30:-1;
                title = "test_model",
                estimator = "ATT"
            )
            
            # Test model properties
            @test model isa CIC
            @test model.title == "test_model"
            @test model.id == :fips
            @test model.t == :day
            @test model.treatment == :gub
            @test model.outcome == :death_rte
            @test model.covariates == matching_covariates
            @test model.timevary == timevary
            @test model.F == 10:20
            @test model.L == -30:-1
            @test model.estimator == "ATT"
            @test length(model.observations) > 0
            @test length(model.ids) > 0
            @test length(model.matches) == length(model.observations)
        end
        
        @testset "Model construction with invalid inputs" begin
            matching_covariates = [:pop_dens, :cumul_death_rate]
            timevary = Dict(:pop_dens => false, :cumul_death_rate => true)
            
            # Test with empty DataFrame - should handle gracefully or error meaningfully
            empty_data = DataFrame()
            @test_throws Exception makemodel(
                empty_data, :day, :fips, :gub, :death_rte,
                matching_covariates, timevary,
                10:20, -30:-1
            )
            
            # Test with missing columns - should error meaningfully
            incomplete_data = test_data[:, [:date, :fips]]  # Missing required columns
            @test_throws Exception makemodel(
                incomplete_data, :day, :fips, :gub, :death_rte,
                matching_covariates, timevary,
                10:20, -30:-1
            )
        end
        
        @testset "Multiple outcome variables" begin
            matching_covariates = [:pop_dens]
            timevary = Dict(:pop_dens => false)
            
            # Test with vector of outcomes
            model_multi = makemodel(
                test_data, :day, :fips, :gub, [:death_rte, :cumul_death_rate],
                matching_covariates, timevary,
                10:20, -30:-1
            )
            
            @test model_multi.outcome == [:death_rte, :cumul_death_rate]
        end
        
        @testset "Edge cases for time ranges" begin
            matching_covariates = [:pop_dens]
            timevary = Dict(:pop_dens => false)
            
            # Test with single time point ranges
            model_single_f = makemodel(
                test_data, :day, :fips, :gub, :death_rte,
                matching_covariates, timevary,
                10:10, -5:-1  # Single F period
            )
            @test length(model_single_f.F) == 1
            
            # Test with overlapping ranges (should be allowed)
            model_overlap = makemodel(
                test_data, :day, :fips, :gub, :death_rte,
                matching_covariates, timevary,
                -5:5, -10:10  # Overlapping F and L
            )
            @test !isempty(intersect(model_overlap.F, model_overlap.L))
        end
    end
    
    @testset "Type system" begin
        @testset "Abstract type hierarchy" begin
            # Test type relationships
            @test CIC <: AbstractCICModel
            @test CICStratified <: AbstractCICModelStratified
            @test AbstractCICModel <: VeryAbstractCICModel
            @test AbstractCICModelStratified <: VeryAbstractCICModel
            
            # Test caliper and refined types
            @test CaliperCIC <: AbstractCICModel
            @test RefinedCIC <: AbstractCICModel
            @test CaliperCICStratified <: AbstractCICModelStratified
            @test RefinedCICStratified <: AbstractCICModelStratified
        end
        
        @testset "TreatmentObservationMatches structure types" begin
            # Test Tob structure creation
            mus = rand(Bool, 10, 5)
            distances = rand(5, 5)
            ranks = Dict(1 => [1, 2, 3], 2 => [2, 1, 3])
            
            tob = TreatmentObservationMatches(mus = mus, distances = distances, ranks = ranks)
            @test tob.eligible_matches == mus
            @test tob.distances == distances
            @test tob.match_rankings == ranks
            
            # Test TreatmentObservationRefinedMatches (Caliper) structure
            tobc = TreatmentObservationRefinedMatches(mus = mus, ranks = ranks)
            @test tobc.eligible_matches == mus
            @test tobc.match_rankings == ranks
            
            # Test TreatmentObservationRefinedMatches (Refined) structure
            tobr = TreatmentObservationRefinedMatches(mus = mus)
            @test tobr.eligible_matches == mus
        end
        
        @testset "Overall results structure" begin
            # Test Overall structure creation
            overall_result = Overall(
                att = 0.5,
                percentiles = [0.1, 0.5, 0.9],
                bayesfactor = 2.0,
                ntreatedmean = 10.0,
                pvalue = 0.05,
                stratum = missing
            )
            
            @test overall_result.att == 0.5
            @test overall_result.percentiles == [0.1, 0.5, 0.9]
            @test overall_result.bayesfactor == 2.0
            @test overall_result.ntreatedmean == 10.0
            @test overall_result.pvalue == 0.05
            @test ismissing(overall_result.stratum)
            
            # Test overall() constructor function
            overall_constructed = overall(
                att = 1.0,
                percentiles = [0.025, 0.975],
                bayesfactor = 3.0
            )
            
            @test overall_constructed.att == 1.0
            @test overall_constructed.percentiles == [0.025, 0.975]
            @test overall_constructed.bayesfactor == 3.0
        end
    end
end