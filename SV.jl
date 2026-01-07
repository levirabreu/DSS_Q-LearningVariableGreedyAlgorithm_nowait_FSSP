# functions to implement the algorithms by Schaller and Valente


# Computing offsets
function compute_offset(N::Int64, M::Int64, P::Matrix{Int64})

    offset = Matrix{Int64}(undef, N, N)
    for j=1:N
        sum_j = Int64[ sum(P[1:i,j]) for i=1:M ]
        for k=1:N
            if j!=k
                sum_l = Int64[ sum(P[1:(i-1),k]) for i=1:M ]
                offset[j,k] = maximum([ sum_j[i] - sum_l[i] for i=1:M ])
            else
                offset[j,k] = 0
            end
        end
    end
    return offset
end


function compute_ET( P::Matrix{Int64},d::Array{Int64}, offset::Matrix{Int64}, S::Array{Int64})

    # computes completion time of each job
    st_1 = Int64[0 for j in eachindex(S)]
    ct_m = Int64[sum(P[1:size(P,1),S[j]]) for j in eachindex(S)]  

    for j in eachindex(S)[2:end]
        st_1[j] = st_1[j-1] + offset[S[j-1],S[j]]
        ct_m[j] += st_1[j] 
    end
   
    # computes ET function
    return sum([max(ct_m[j] - d[job],0) for (j,job) in enumerate(S)]) + sum([ max(d[job] - ct_m[j],0) for (j,job) in enumerate(S)])

end



function LIN1(N::Int64,M::Int64,P::Matrix{Int64},d::Array{Int64}, w::Float64)

    S = Int64[]  # scheduled
    U = Int64[j for j=1:N]   # unscheduled
    
    st_1 = 0    # starting time of the last jobs in the first machine
    ct_m = 0   # completion time of the last job 
    
    offset = compute_offset(N, M, P)
    
    while length(U) > 1
    
        indicator = [0.0 for j=1:length(U)]
        # calculating threshold
        lowest_offset = maximum(offset[:,:])    # calculate lowest offset of unscheduled
        for j in U, k in U
            if k!=j && offset[j,k] < lowest_offset
                lowest_offset = offset[j,k]
            end
        end
        if length(S) == 0
            slk_thr = w * (sum(P[1,j] for j in U) + (N-length(S)-1)* lowest_offset )
        else
            slk_thr = w * (sum(P[1,j] for j in U) + (N-length(S)-1)* lowest_offset + minimum([offset[S[length(S)],j] for j in U]) )
        end
    
        for (j,job) in enumerate(U) # for each unscheduled job
    
            # compute P_j and slacks for job j
            if length(S) == 0
                st_1_j = 0
                ct_m_j = sum(P[1:M,job])
                P_j = sum(P[1:M,job])    
            else
                st_1_j = st_1 + offset[S[length(S)],job]
                ct_m_j = st_1_j + sum(P[1:M,job])
                P_j = offset[S[length(S)],job] + sum(P[1:M,job])
            end
            slack_j = d[job] - ct_m_j
    
            # compute the indicator
            if slack_j <= 0
                indicator[j] = 1/P_j
            elseif slack_j > 0 && slk_thr > slack_j
                indicator[j] = 1/P_j - slack_j/slk_thr * 2/P_j
            else
                indicator[j] = (-1)*1/P_j
            end
    
        end
    
        selected_index = sortperm(indicator, rev=true)[1]
   
        # update starting times and completion times
        if length(S) == 0
            st_1 = 0
            ct_m = sum(P[1:M,U[selected_index]])
        else
            st_1 = st_1 + offset[S[length(S)],U[selected_index]]
            ct_m = st_1 + sum(P[1:M,U[selected_index]])
        end
        # update scheduled and unscheduled jobs
        push!(S, U[selected_index])
        deleteat!(U, selected_index)
        
    end

    push!(S,U[1])   # last job
    
    return S

end


# H1 is a modification of LIN1 based on two different ways to compute P_j (previosly computed as P_j = offset[S[length(S)],j] + sum(P[1:M,j]))
# 1) only the processing time in the first machine is considered in the second term in P_j
# 2) A third term based on the average offset of the non scheduled jobs is added to P_j
# As a result of these changes:
# P_j = offset[S[length(S)],j] + P[1,j] + AD_j if P_j >0, else P_j = 1
# where AD_j = offset[S[length(S)],j] - sum( offset[S[length(S)],j] for k in U)/(N - length(S)-1) 
function H1(N::Int64,M::Int64,P::Matrix{Int64},d::Array{Int64}, w::Float64)

    S = Int64[]  # scheduled
    U = Int64[j for j=1:N]   # unscheduled
    
    st_1 = 0    # starting time of the last jobs in the first machine
    ct_m = 0   # completion time of the last job 
    
    offset = compute_offset(N, M, P)
    
    while length(U) > 1
    
        indicator = Float64[0.0 for j=1:length(U)]
        # calculating threshold
        lowest_offset = maximum(offset[:,:])    # calculate lowest offset of unscheduled
        for j in U
            for k in U
                if k!=j && offset[j,k] < lowest_offset
                    lowest_offset = offset[j,k]
                end
            end
        end
        if length(S) == 0
            slk_thr = w * (sum(P[1,j] for j in U) + (N-length(S)-1)* lowest_offset )
        else
            slk_thr = w * (sum(P[1,j] for j in U) + (N-length(S)-1)* lowest_offset + minimum([offset[S[length(S)],j] for j in U]) )
        end
    
        for (j,job) in enumerate(U) # for each unscheduled job
    
            # compute P_j and slacks for job j
            if length(S) == 0
                st_1_j = 0
                ct_m_j = sum(P[1:M,job])
                P_j = P[1,job]
            else
                st_1_j = st_1 + offset[S[length(S)],job]
                ct_m_j = st_1_j + sum(P[1:M,job])
                AD_j = offset[S[length(S)],j] - sum( offset[S[length(S)],j] for k in U)/(N - length(S)-1)
                P_j = offset[S[length(S)],j] + P[1,j] + AD_j
                if P_j < 1
                    P_j = 1
                end
            end
            slack_j = d[job] - ct_m_j
    
            # compute the indicator
            if slack_j <= 0
                indicator[j] = 1/P_j
            elseif slack_j > 0 && slk_thr > slack_j
                indicator[j] = 1/P_j - slack_j/slk_thr * 2/P_j
            else
                indicator[j] = (-1)*1/P_j
            end
    
        end
    
        selected_index = sortperm(indicator, rev=true)[1]
   
        # update starting times and completion times
        if length(S) == 0
            st_1 = 0
            ct_m = sum(P[1:M,U[selected_index]])
        else
            st_1 = st_1 + offset[S[length(S)],U[selected_index]]
            ct_m = st_1 + sum(P[1:M,U[selected_index]])
        end
        # update scheduled and unscheduled jobs
        push!(S, U[selected_index])
        deleteat!(U, selected_index)
        
    end

    push!(S,U[1])   # last job
    
    return S

end

# Takes the H1 heuristic and applies a local search to improve the solution
# The local search consists of removing a job and inserting it in another position
function H1INS(N::Int64,M::Int64,P::Matrix{Int64},d::Array{Int64}, w::Float64)

    offset = compute_offset(N, M, P)

    # initial solution
    incumbent_S = H1(N,M,P,d, w)
    incumbent_ET = compute_ET(P,d, offset, incumbent_S)
    #print("\nInitial solution H1: $incumbent_ET - $incumbent_S\n") 

    improvement = true
    n_iterations = 0
    while improvement 
        improvement = false
        for j=1:N, k=1:N
            if j != k
                # insert job j in position k    
                S = copy(incumbent_S)
                to_insert = S[j]
                deleteat!(S, j)
                insert!(S, k, to_insert)
                ET = compute_ET(P,d, offset, S)
                if ET < incumbent_ET
                    #print("Iteration $n_iterations: $incumbent_ET -> $ET\n")
                    incumbent_S = copy(S)
                    incumbent_ET = ET
                    improvement = true
                    break
                end
            end
        end
        n_iterations += 1
    end

    return incumbent_S

end 