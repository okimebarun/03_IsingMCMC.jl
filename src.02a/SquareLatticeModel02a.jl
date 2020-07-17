module SquareLatticeModel
##export ParamDict, getNeighbors, genPoly

mutable struct ParamDict
    nrows::Int
    ncols::Int
    N::Int
end

function genFunctor(d)
    d = d
    
    function roundN(n, N)
        m = (n - 1 + N) % N
        m + 1
    end

    function getXY(k, d=d)
        x = roundN(k, d.ncols)
        y = roundN(div(k-1, d.ncols)+1, d.nrows)
        (x, y)
    end

    function getK( x, y, d=d)
        # must round each x, y
        xn = (x - 1 + d.ncols) % d.ncols + 1
        yn = (y - 1 + d.nrows) % d.nrows + 1
        #
        k = (yn - 1) * d.ncols + xn
    end

    function getNeighbors( k, d=d)
        (x, y) = getXY(k,d)
        [getK(x-1, y, d), getK(x, y-1, d), getK(x+1, y, d), getK(x, y+1, d)]
    end

    # polygon for each cell
    function genPoly(k, d)
        (x, y) = getXY(k, d)
        w = 0.5
        h = 0.5
        (x .+ [-w,  w, w, -w, -w],
         y .+ [-h, -h, h,  h, -h])
    end

    return (getNeighbors, genPoly)
end # of functor

end # of module