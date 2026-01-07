function calculate_temperature(P,tal)
    m = size(P,1)
    n = size(P,2)
    T = tal*(sum(P)/(10*m*n))
    return(T)
end
