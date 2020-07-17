module HexaLatticeModel
##export ParamDict, getNeighbors, genPoly

mutable struct ParamDict
    nrows::Int
    ncols::Int
    N::Int
end

function genFunctor(d)
    d = d
    # cell size
    w = 0.5
    h = 2w/sqrt(3)

    function roundN(n, N)
        m = (n - 1 + N) % N
        m + 1
    end

    function getXY(k, d=d)
        m = roundN(k, d.ncols)
        n = roundN(div(k-1, d.ncols)+1, d.nrows)
        if n % 2 == 1
            x = m + 0.5
            y = n*(1.5h)
        else
            x = m
            y = n*(1.5h)
        end
        (x, y)
    end

    function getK( x, y, d=d)
        # must round each x, y
        n = (round(Int,y/(1.5h)) - 1 + d.nrows) % d.nrows + 1
        #
        if n % 2 == 1
            xn = Int(x - 0.5)
        else
            xn = Int(x)
        end
        m = (xn - 1 + d.ncols) % d.ncols + 1
        k = (n - 1) * d.ncols + m
    end

    function getNeighbors( k, d=d)
        (x, y) = getXY(k, d)
        [getK(x-2w,      y, d), getK(x+2w,      y, d),
         getK(x- w, y-1.5h, d), getK(x+ w, y-1.5h, d),
         getK(x- w, y+1.5h, d), getK(x+ w, y+1.5h, d)]
    end

    # polygon for each cell
    function genPoly(k, d=d)
        (x, y) = getXY(k, d)
        (x .+ [  -w, 0,   w,    w,  0,   -w,  -w],
         y .+ [ h/2, h, h/2, -h/2, -h, -h/2, h/2])
    end
    
    return (getNeighbors, genPoly)
end # of functor

end # of module