using Random
function random_initial_solution(n)
    random_pi = zeros(Int64,n)
    random_pi=randperm(n)
    return(random_pi)
end