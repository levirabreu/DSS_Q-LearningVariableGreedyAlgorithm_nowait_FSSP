function earliness_tardiness(Pi,P,d, M, N)
    m = size(P,1)
    n = length(Pi)
    N = size(P,2)
    dist = zeros(Int64,m)
    C = zeros(Int64,N)
    T = zeros(Int64,N)
    D = zeros(Int64,n)
    aux1 = zeros(Int64,N)
    aux2 = zeros(Int64,N)
    P_ = zeros(Int64,N)

    global s1 = 0
    global s2 = 0

    j = 2
    k = 1

    for a in 1:n-1
        for i in 1:m
            for h in 1:i
                s1 = s1 + P[h,Pi[a]]
                #println(h,"  ",s1)
            end
            if (i > 1)
                for h in 1:i-1
                    s2 = s2 + P[h,Pi[a+1]]
                end
            else
                s2 = 0
            end
            dist[i] = s1 - s2
            D[a] = maximum(dist)
            s1 = 0
            s2 = 0
        end
    end

    for j in 1:n
        P_[j]=sum(P[:,j])
    end

    for j in 2:n
        aux1[j] += D[j-1]
    end

    aux2=cumsum(aux1)

    for j in 1:n
        C[Pi[j]]= P_[Pi[j]] + aux2[j]
    end

    for j in 1:N
        if true#(C[j]-d[j]) > 0
            T[j] = abs(C[j]-d[j])
        end
    end
    
    #println(C)

    return(sum(T))
end

 #p, d = read_my_file("C:/Users/levi1/OneDrive - Universidade Federal do Ceará/UFC/Projetos UFC/flowshop/Flow shop no-wait e ET/Test instances/FSPNWDD/1-100 mgt3/19.txt")

 #M = 5
 #N = 20

 #sol = [31, 28, 18, 21, 25, 34, 36, 17, 11, 14, 37, 13, 7, 26, 32, 33, 20, 15, 12, 22, 16, 27, 6, 8, 38, 24, 29, 23, 9, 1, 4, 10, 40, 3, 35, 2, 39, 30, 19, 5]

"""
p = [[2, 2, 3],
              [3, 4, 2],
              [2, 3, 4]]
p = mapreduce(permutedims, vcat, p)

d = [8, 10, 12]

m, n = 3, 3

sol = [1, 3, 2]
"""

#ET = earliness_tardiness(sol,p,d)

