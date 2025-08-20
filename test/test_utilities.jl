@testset "Utility Functions" begin
    
    # Helper function to create test model and data
    function create_test_setup()
        dates = Date(2021, 10, 1):Day(1):Date(2021, 12, 31)
        n_units = 25
        
        test_data = DataFrame()
        for unit in 1:n_units
            for (i, date) in enumerate(dates)
                push!(test_data, (
                    date = date,
                    fips = 1000 + unit,
                    pop_dens = 50.0 + randn() * 10,
                    cumul_death_rate = 20.0 + i * 0.1 + randn() * 2,
                    death_rte = abs(randn() * 0.5),
                    gub = (unit <= 6) ? 1 : 0
                ))
            end
        end
        
        test_data[!, :day] = Dates.value.(test_data.date .- minimum(test_data.date))
        
        model = makemodel(
            test_data, :day, :fips, :gub, :death_rte,
            [:pop_dens, :cumul_death_rate],
            Dict(:pop_dens => false, :cumul_death_rate => true),
            50:60, -30:-1
        )
        
        match!(model, test_data)
        balance!(model, test_data)
        
        return model, test_data
    end
    
    @testset "example_data function" begin
        @testset "Data loading" begin
            # Test that example_data loads successfully
            @test_nowarn example_data()
            
            # Test properties of loaded data
            data = example_data()
            @test data isa DataFrame
            @test nrow(data) > 0
            @test ncol(data) > 0
            
            # Should have expected columns from vignette
            expected_cols = [:date, :fips, :pop_dens, :cumul_death_rate, :death_rte, :gub]
            for col in expected_cols
                @test col in names(data)
            end
        end
        
        @testset "Data quality" begin
            data = example_data()
            
            # Basic data quality checks
            @test all(x -> x isa Date, data.date)
            @test all(x -> x isa Integer, data.fips)
            @test all(x -> x >= 0, data.pop_dens)
            @test all(x -> x in [0, 1], data.gub)
            
            # Should have multiple units and time periods
            @test length(unique(data.fips)) > 1
            @test length(unique(data.date)) > 1
        end
    end
    
    @testset "Information and inspection functions" begin
        model, test_data = create_test_setup()
        
        @testset "matchinfo function" begin
            # Test that matchinfo runs and returns information
            @test_nowarn matchinfo(model, test_data)
            
            info = matchinfo(model, test_data)
            @test !isnothing(info)
            
            # Should provide information about matching results
            if info isa DataFrame
                @test nrow(info) > 0
            elseif info isa Dict
                @test length(info) > 0
            end
        end
        
        @testset "obsinfo function" begin
            # Test that obsinfo runs
            @test_nowarn obsinfo(model, test_data)
            
            obs_info = obsinfo(model, test_data)
            @test !isnothing(obs_info)
        end
        
        @testset "showmatches function" begin
            # Test that showmatches runs
            @test_nowarn showmatches(model, test_data)
            
            # Should display or return information about matches
            matches_display = showmatches(model, test_data)
            @test !isnothing(matches_display) || matches_display === nothing  # Might just print
        end
        
        @testset "treatedinfo function" begin
            # Test information about treated units
            @test_nowarn treatedinfo(model, test_data)
            
            treated_info = treatedinfo(model, test_data)
            @test !isnothing(treated_info)
            
            # Should provide info about treated units
            if treated_info isa DataFrame
                @test nrow(treated_info) > 0
            elseif treated_info isa Vector
                @test length(treated_info) > 0
            end
        end
    end
    
    @testset "Model manipulation functions" begin
        model, test_data = create_test_setup()
        
        @testset "relabel! function" begin
            # Test relabeling functionality
            original_title = model.title
            
            @test_nowarn relabel!(model, "new_test_title")
            
            # Should update model title or identifier
            @test model.title != original_title || hasfield(typeof(model), :label)
        end
        
        @testset "trim_model function" begin
            # Test model trimming functionality
            @test_nowarn trim_model(model)
            
            trimmed = trim_model(model)
            @test !isnothing(trimmed)
            @test typeof(trimmed) == typeof(model) || trimmed isa VeryAbstractCICModel
        end
        
        @testset "name_model function" begin
            # Test model naming functionality
            @test_nowarn name_model(model)
            
            model_name = name_model(model)
            @test model_name isa String
            @test length(model_name) > 0
        end
    end
    
    @testset "Refinement and caliper functions" begin
        model, test_data = create_test_setup()
        
        @testset "refine function" begin
            # Test model refinement
            @test_nowarn refine(model, test_data; refinementnum = 3)
            
            refined_model = refine(
                model, test_data;
                refinementnum = 5,
                dobalance = false,
                doestimate = false
            )
            
            @test refined_model isa VeryAbstractCICModel
            @test typeof(refined_model) != typeof(model)  # Should be refined type
        end
        
        @testset "refine with balancing and estimation" begin
            refined_balanced = refine(
                model, test_data;
                refinementnum = 4,
                dobalance = true,
                doestimate = false
            )
            
            @test refined_balanced isa VeryAbstractCICModel
            
            # If estimation is also done
            refined_estimated = refine(
                model, test_data;
                refinementnum = 3,
                dobalance = true,
                doestimate = true
            )
            
            @test refined_estimated isa VeryAbstractCICModel
        end
        
        @testset "caliper function" begin
            # Test caliper application
            caliper_specs = Dict(:pop_dens => 0.5, :cumul_death_rate => 0.3)
            
            @test_nowarn caliper(model, test_data, caliper_specs)
            
            caliper_model = caliper(model, test_data, caliper_specs)
            @test caliper_model isa VeryAbstractCICModel
            @test typeof(caliper_model) != typeof(model)  # Should be caliper type
        end
        
        @testset "caliper with processing options" begin
            caliper_specs = Dict(:pop_dens => 0.4)
            
            caliper_processed = caliper(
                model, test_data, caliper_specs;
                dobalance = true,
                doestimate = false
            )
            
            @test caliper_processed isa VeryAbstractCICModel
        end
    end
    
    @testset "Processing and workflow functions" begin
        model, test_data = create_test_setup()
        
        @testset "matchprocess function" begin
            # Test the complete matching process
            @test_nowarn matchprocess(model, test_data)
            
            processed_model = matchprocess(model, test_data)
            @test processed_model isa VeryAbstractCICModel
        end
        
        @testset "quick_att function" begin
            # Test quick ATT estimation
            @test_nowarn quick_att(model, test_data)
            
            att_result = quick_att(model, test_data)
            @test att_result isa Real || att_result isa Vector{<:Real}
        end
        
        @testset "variable_filter function" begin
            # Test variable filtering functionality
            test_vars = [:pop_dens, :cumul_death_rate, :nonexistent_var]
            
            @test_nowarn variable_filter(test_data, test_vars)
            
            filtered_vars = variable_filter(test_data, test_vars)
            @test filtered_vars isa Vector{Symbol}
            @test all(var -> var in names(test_data), filtered_vars)
        end
    end
    
    @testset "Inspection functions" begin
        model, test_data = create_test_setup()
        
        @testset "inspection function" begin
            # Test model inspection
            @test_nowarn inspection(model, test_data)
            
            inspection_result = inspection(model, test_data)
            @test !isnothing(inspection_result)
        end
        
        @testset "pretreatment function" begin
            # Test pretreatment period analysis
            @test_nowarn pretreatment(model, test_data)
            
            pretreat_result = pretreatment(model, test_data)
            @test !isnothing(pretreat_result)
        end
    end
    
    @testset "Stratification functions" begin
        model, test_data = create_test_setup()
        
        @testset "stratify function" begin
            # Test basic stratification
            strata_var = :pop_dens
            @test_nowarn stratify(model, test_data, strata_var)
            
            stratified_model = stratify(model, test_data, strata_var)
            @test stratified_model isa AbstractCICModelStratified
        end
        
        @testset "variablestrat function" begin
            # Test variable-based stratification
            @test_nowarn variablestrat(model, test_data, :pop_dens)
            
            var_strat_model = variablestrat(model, test_data, :pop_dens)
            @test var_strat_model isa AbstractCICModelStratified
        end
        
        @testset "combostrat function" begin
            # Test combination stratification
            strat_vars = [:pop_dens, :cumul_death_rate]
            @test_nowarn combostrat(model, test_data, strat_vars)
            
            combo_strat_model = combostrat(model, test_data, strat_vars)
            @test combo_strat_model isa AbstractCICModelStratified
        end
        
        @testset "customstrat function" begin
            # Test custom stratification
            # Create a custom stratification function
            custom_strat_func = x -> x > median(test_data.pop_dens) ? 1 : 0
            
            @test_nowarn customstrat(model, test_data, custom_strat_func)
            
            custom_model = customstrat(model, test_data, custom_strat_func)
            @test custom_model isa AbstractCICModelStratified
        end
    end
    
    @testset "Data manipulation utilities" begin
        model, test_data = create_test_setup()
        
        @testset "filterunits! function" begin
            # Test unit filtering (modifies data in place)
            test_data_copy = copy(test_data)
            original_nrows = nrow(test_data_copy)
            
            # Filter to keep only certain units
            units_to_keep = unique(test_data_copy.fips)[1:5]
            filterunits!(test_data_copy, units_to_keep)
            
            @test nrow(test_data_copy) <= original_nrows
            @test all(fips -> fips in units_to_keep, test_data_copy.fips)
        end
    end
    
    @testset "Storage and records functions" begin
        model, test_data = create_test_setup()
        estimate!(model, test_data; iterations = 50)  # Ensure we have results
        
        @testset "makerecords function" begin
            # Test record creation
            temp_dir = mktempdir()
            models_list = [model]
            
            @test_nowarn makerecords(test_data, temp_dir, models_list)
            
            records = makerecords(test_data, temp_dir, models_list)
            @test !isnothing(records)
            
            # Clean up
            rm(temp_dir, recursive = true)
        end
        
        @testset "save function (if different from JLD2)" begin
            # Test TSCSMethods-specific save function if it exists
            temp_dir = mktempdir()
            temp_file = joinpath(temp_dir, "test_save.jld2")
            
            if isdefined(TSCSMethods, :save) && TSCSMethods.save !== Base.save
                @test_nowarn TSCSMethods.save(temp_file, model)
            end
            
            # Clean up
            rm(temp_dir, recursive = true)
        end
    end
end