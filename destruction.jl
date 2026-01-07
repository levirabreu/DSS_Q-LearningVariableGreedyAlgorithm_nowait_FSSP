using StatsBase
function destruction(pi_0,k_destruction)
    pi_D = sample(1:length(pi_0),k_destruction, replace = false)
    pi_P = setdiff(pi_0, pi_D)
    return(pi_D, pi_P)
end
