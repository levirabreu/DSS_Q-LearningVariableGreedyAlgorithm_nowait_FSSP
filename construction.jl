function construction(pi_D, pi_P,M,N,P,d)
    best_Pi = []
    global best_earliness_tardiness = 100000000000
    current_earliness_tardiness = 0
    best_position = 0
    best_value = 0
    for k in 1:lastindex(pi_D)
        best_earliness_tardiness = 100000000000
        for l in 1:lastindex(pi_P)+1
            insert!(pi_P,l,pi_D[k])
            current_earliness_tardiness = earliness_tardiness(pi_P,P,d,M,k)
            if current_earliness_tardiness < best_earliness_tardiness
                best_earliness_tardiness = current_earliness_tardiness
                best_Pi = copy(pi_P)
                best_position = l
                best_value = pi_D[k]
                deleteat!(pi_P,l)
            else
                deleteat!(pi_P,l)
            end
        end
        insert!(pi_P,best_position, best_value)
    end

    best_earliness_tardiness = earliness_tardiness(pi_P,P,d,M,N)
    return(pi_P,best_earliness_tardiness)
end
