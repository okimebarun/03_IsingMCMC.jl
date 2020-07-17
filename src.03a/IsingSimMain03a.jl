# load packages
using Plots
#pyplot()
gr()

# load modules
# select a model later
include("HexaLatticeModel03a.jl")
include("SquareLatticeModel03a.jl")
# MCMC
include("IsingMCMC03a.jl")

###################################################################################################
# view functions
function drawEnergy(simE, trial, T, d=d, N=N)
    xs = [t for t in 0:1000:trial]
    ys = simE
    p = Plots.plot(xs, ys, label="Energy history with T=$(T)")
end

function drawSim(Sf, d=d, N=N)
    snap = Sf
    
    p=Plots.plot([],[],legend=false, size=(500,500))
    
    npoly = length(genPoly(1,d)[1])
    nbuff = (npoly+1)*N
    p0xs = fill(NaN, nbuff);p0ys = fill(NaN, nbuff)
    p1xs = fill(NaN, nbuff);p1ys = fill(NaN, nbuff)
    i1 = 1; i0 = 1
    @inbounds @simd for i in 1:N
        # bits
        (xs,ys) = genPoly(i, d)
        if (snap[i] == 1) # up-spin
            p1xs[i1:(i1+npoly-1)] .= xs
            p1ys[i1:(i1+npoly-1)] .= ys
            i1 += npoly+1
        else # down-spin
            p0xs[i0:(i0+npoly-1)] .= xs
            p0ys[i0:(i0+npoly-1)] .= ys
            i0 += npoly+1
        end
    end    
    poly0 = Shape(p0xs, p0ys) # creating polygon
    p = Plots.plot!(poly0, fillcolor = plot_color(:blue)  ,legend=false)
    poly1 = Shape(p1xs, p1ys)
    p = Plots.plot!(poly1, fillcolor = plot_color(:yellow),legend=false)
end

###################################################################################################
# prepare
# setting for energy function

function genInitialJB(d=d)
    Jd = [Dict{Int,Int8}() for i in 1:d.N]
    Bd = zeros(Int8, d.N)
    (Jd, Bd)
end

function setJB!( Jd=Jd, Bd=Bd, d=d)
    J0 = 1
    B0 = 0
    for i in 1:d.N
        for j in getNeighbors(i,d)
            @inbounds Jd[i][j] = J0
        end
        @inbounds Bd[i] = B0
    end
end

###################################################################################################
# do simulation

# setting for model
d = SquareLatticeModel.ParamDict(128,128,0)
#d = SquareLatticeModel.ParamDict(64,64,0)
#d = SquareLatticeModel.ParamDict(32,32,0)
#d = SquareLatticeModel.ParamDict(8,8,0)
N = d.nrows * d.ncols
d.N = N

# main function
# mdlno: 1(square), 2(hexa)
# T: temperature e.g. T=2.26
function main(mdlno, T) 
    # setting for simulation
    trial = 100000 # MCMC trial
    
    # gen functors as global
    global getNeighbors, genPoly
    if mdlno == 1
        (getNeighbors, genPoly)= SquareLatticeModel.genFunctor(d)
    elseif mdlno == 2
        (getNeighbors, genPoly)= HexaLatticeModel.genFunctor(d)
    end
    
    # generate Jd, Bd
    (Jd, Bd) = genInitialJB(d)
    setJB!(Jd, Bd)
    
    # simulation
    global E, MCMC
    (E, MCMC) = IsingMCMC.genFunctor(Jd, Bd)
    print("start MCMC...")
    @time (simE, Sf) = MCMC(T, N, trial);

    # save images
    # energy history
    print("start drawEnergy...")
    @time p = drawEnergy(simE, trial, T)
    Plots.savefig(p, "Ising_energy_M$(mdlno)_T$(T).png")

    # cell state
    print("start drawSim...")
    @time p = drawSim(Sf, d)
    Plots.savefig(p, "Ising_sim_M$(mdlno)_T$(T).png")
    print("done.")
end
###################################################################################################
