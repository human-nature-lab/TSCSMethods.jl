@testset "Simple Validation Test" begin
    # Test very basic functionality
    
    @testset "Basic type validation" begin
        # Test that makemodel requires DataFrame
        @test_throws MethodError makemodel("not_a_dataframe", :t, :id, :treat, :out, [:cov], Dict(:cov => false), 1:5, -5:-1)
        
        # Test that empty DataFrame is caught
        @test_throws ArgumentError makemodel(DataFrame(), :t, :id, :treat, :out, [:cov], Dict(:cov => false), 1:5, -5:-1)
        
        println("Basic validation tests passed!")
    end
    
    @testset "Simple working example" begin
        # Create minimal working data
        test_data = DataFrame(
            time_var = [1, 2, 3, 4, 5, 6],
            unit_id = [1, 1, 1, 2, 2, 2],
            treatment = [0, 0, 1, 0, 0, 0],
            outcome = [1.0, 1.1, 1.2, 2.0, 2.1, 2.2],
            covariate1 = [10.0, 10.1, 10.2, 20.0, 20.1, 20.2]
        )
        
        println("Test data columns: ", names(test_data))
        
        # This should work
        @test_nowarn model = makemodel(
            test_data, :time_var, :unit_id, :treatment, :outcome,
            [:covariate1], Dict(:covariate1 => false),
            1:2, -2:-1
        )
        
        println("Simple working example passed!")
    end
end