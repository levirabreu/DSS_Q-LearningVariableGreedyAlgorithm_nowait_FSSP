### Test QLVG ###

include("calculate_temperature.jl")
include("earliness_tardiness.jl")
include("destruction.jl")
include("construction.jl")
include("random_initial_solution.jl")
include("QLVG.jl")
include("SV.jl")

P = [ 1 68 9 30 70 33 49 98 15 1 56 65 9 58 71 15 69 88 18 50;
      4 32 48 92 85 47 15 50 15 78 21 100 78 96 99 92 29 17 16 29;
     87 17 8 37 72 25 88 89 51 66 69 25 50 69 78 79 10 28 27 16;
     21 38 85 78 31 83 29 83 3 78 60 68 89 34 75 59 63 55 18 61;
     28 43 6 33 17 28 78 3 60 71 96 30 51 1 20 87 29 97 58 74]

d = [1038, 905, 905, 901, 911, 1048, 1050, 965, 1011, 976, 1037, 988, 957, 1052, 1052, 918, 908, 903, 901, 928]

N = length(d)
M = size(P,1) 
time_limit = N*(M/2)*60/1000

size_prop = 0.7

best_pi, best_jit = run_QLVG(P, d, size_prop, M, N, time_limit)