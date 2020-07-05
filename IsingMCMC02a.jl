module IsingMCMC
##export E, MCMC

function genFunctor( getJ, getB)
    getJ = getJ
    getB = getB
    
    # E(S) = - 1/2\sum_{i,j, i not j} J_{i,j} si sj - \sum_i B_i si
    function E(S, N)
        Et = 0
        for i in 1:N
            J = Main.getNeighbors(i)
            Et += -0.5*S[i]*sum( [getJ(i,j) * S[j] for j in J] )
            Et += -    S[i]*getB(i)
        end
        Et
    end

    # differential energy
    function dE(S, k, N)
        dEt = 0
        J = Main.getNeighbors(k)
        dEt +=   S[k]*sum( [(getJ(k,j) + getJ(j,k)) * S[j] for j in J] )
        dEt += 2*S[k]*getB(k)
        dEt
    end

    # probability for dE
    PrE(dE, T) = exp(-dE/T)

    # flip spin at x
    function flipx!(list, x)
        list[x] *= -1
    end

    # MCMC
    function MCMC( T, N, trial)
        # initialize
        sim = []
        S = ones(N)
        for i in 1:N
            if rand() < 0.5
                flipx!(S,i) # random flip at first
            end
        end
        # Gibbs sampling position
        gpos = [rand(1:N) for i in 1:trial]
        # random value in [0,1]
        rn = [rand() for i in 1:trial]

        # MCMC trial
        for t in 1:trial
            k = gpos[t]
            de = dE(S, k, N)
            if rn[t] < PrE(de, T) # MH criteria 
                flipx!(S, k) # change k
            end
            push!(sim, copy(S))
        end

        sim
    end
    
    return (E, MCMC)
end # of functor

end # of module