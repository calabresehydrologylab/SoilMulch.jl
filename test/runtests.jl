using SoilMulch
using Test



@testset "SoilMulch.jl" begin
    # Write your tests here.
    # plant_transpiration(1, 0.15, 0.5, 1, 0.5)
    @test plant_transpiration(1, 0.15, 0.5, 1, 0.5) == 0.5
    # @test soil_leakage(1, 150, 14) == 150
    @test soil_leakage(1, 250, 14) == 250
    @test soil_evaporation(1, 0.15, 1, 0.5) == 0.5
end
