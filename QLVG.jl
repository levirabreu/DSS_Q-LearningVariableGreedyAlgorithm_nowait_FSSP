function run_QLVG(P, d, size_prop, M, N, time_limit)
    elapsed = 0
    t0 = time()
    current_pi = H1INS(N,M,P,d, 0.6)
    current_ET = earliness_tardiness(current_pi,P,d, M, N)
    best_ET = current_ET
    best_pi = copy(current_pi)

    A_k = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    Q = Dict(0.1=>0.0, 0.2=>0.0, 0.3=>0.0, 0.4=>0.0, 0.5=>0.0, 0.6=>0.0, 0.7=>0.0, 0.8=>0.0, 0.9=>0.0)

    ϵ = 0.8
    α = 0.6
    γ = 0.8

    size_prop = rand(A_k)

    k = 1
    k_max =  ceil(Int64,size_prop*(length(d)-1))

    while (elapsed <= time_limit)
        pi_D,pi_P = destruction(current_pi, k)
        new_pi, new_ET = construction(pi_D, pi_P,M,N,P,d)

        if new_ET < current_ET
            current_ET = new_ET
            current_pi = copy(new_pi)
            k = 1
            r = 1/new_ET
            Q[size_prop] = Q[size_prop] + α*(r + γ*maximum(Q[a] for a in A_k) - Q[size_prop])
        else
            k += 1
            if (k > k_max)
                k = 1
                current_pi = random_initial_solution(length(d))
                pi_D,pi_P = destruction(current_pi, k_max)
                current_pi, current_ET = construction(pi_D, pi_P,M,N,P,d)
            end
        end

        if current_ET < best_ET
            best_ET = current_ET
            best_pi = copy(current_pi)
        end

        if rand() >= ϵ
            better_idx = argmax(Q[a] for a in A_k)
            size_prop = A_k[better_idx]
            k_max =  ceil(Int64,size_prop*(length(d)-1))
        else
            size_prop = rand(A_k)
            k_max =  ceil(Int64,size_prop*(length(d)-1))
        end


        t1 = time()
        elapsed = t1-t0


        
    end

    return(best_pi, best_ET)
end