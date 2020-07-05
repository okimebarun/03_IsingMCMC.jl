# load packages
using Plots
#pyplot()
gr()

# load modules
# select a model later
include("HexaLatticeModel02a.jl")
include("SquareLatticeModel02a.jl")
# MCMC
include("IsingMCMC02a.jl")

###################################################################################################
# view functions
function drawEnergy(sim, trial, T, d=d, N=N)
    xs = [t for t in 1:1000:trial]
    ys = [E(sim[t],N) for t in xs]
    p = Plots.plot(xs, ys, label="Energy history with T=$(T)")
end

function drawSim(sim, f, d=d, N=N)
    snap = sim[f]
    
    p=Plots.plot([],[],legend=false, size=(500,500))
    
    p0xs = [];p0ys = []
    p1xs = [];p1ys = []
    for i in 1:N
        # bits
        (xs,ys) = genPoly(i, d)
        if (snap[i] == 1) # up-spin
            p1xs = vcat(p1xs, NaN, xs)
            p1ys = vcat(p1ys, NaN, ys)
        else # down-spin
            p0xs = vcat(p0xs, NaN, xs)
            p0ys = vcat(p0ys, NaN, ys)
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
function getJ(i, j, d=d)
    J0 = 1
    (j in getNeighbors(i,d)) ? J0 : 0
end

function getB(i, d=d)
    B0 = 0
    B0
end

###################################################################################################
# do simulation

# setting for model
d = SquareLatticeModel.ParamDict(64,64,0)
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

    # simulation
    global E, MCMC
    (E, MCMC) = IsingMCMC.genFunctor(getJ, getB)
    print("start MCMC...")
    @time sim = MCMC(T, N, trial);

    # save images
    # energy history
    print("start drawEnergy...")
    @time p = drawEnergy(sim, trial, T)
    Plots.savefig(p, "Ising_energy_M$(mdlno)_T$(T).png")

    # cell state
    print("start drawSim...")
    @time p = drawSim(sim, trial, d)
    Plots.savefig(p, "Ising_sim_M$(mdlno)_T$(T).png")
    print("done.")
end
###################################################################################################
