#generate inital north boundary mpses of bond dimension chi
function gen_north_bounds(peps,chi)
    T = eltype(peps[1,1])

    boundaries = map(1:size(peps,1)) do i
        lds = [1]
        rds = [1]

        for j = 1:size(peps,2)
            push!(lds,min(size(peps[i,j],North)*size(peps[i,j],North)*lds[j],chi))
            push!(rds,min(size(peps[i,end-j+1],North)*size(peps[i,end-j+1],North)*rds[j],chi))
        end
        rds = reverse(rds)

        cb = map(1:size(peps,2)) do j
            psize = size(peps[i,j],North)

            if i>1
                lD = min(lds[j],rds[j])
                rD = min(lds[j+1],rds[j+1])
            else
                lD = 1
                rD = 1
            end

            randn(T,lD,psize,psize,rD)
        end
    end
end

function gen_boundaries(peps,chi)
    T = eltype(peps[1,1]);

    bs = map([Dirs...]) do dir
        gen_north_bounds(rotate_north(peps,dir),chi)
    end
end

function calc_north_boundaries(peps,pboundaries;maxiter=Defaults.maxiter,tol=Defaults.tol)
    #use vomps to optimize them
    boundaries = Zygote.bufferfrom(pboundaries)
    for i in 1:(size(peps,1)-1)
        #=
        I want to take a peps - slice but can't (wouldn't be zygote buffer)
        pepsline = peps[i,:]
        =#

        T = typeof(peps[1,1])
        pepsline = Zygote.Buffer(T[])
        for t = 1:size(peps,2)
            push!(pepsline,peps[i,t])
        end

        boundaries[i+1] = north_vomps(boundaries[i],pepsline,boundaries[i+1],tol=tol,maxiter=maxiter)

    end
    return copy(boundaries)
end

function calc_boundaries(peps,old_boundaries;maxiter=Defaults.maxiter,tol=Defaults.tol)
    boundaries = Zygote.Buffer(old_boundaries)

    for cd in Dirs
        boundaries[cd] = calc_north_boundaries(rotate_north(peps,cd),old_boundaries[cd],maxiter=maxiter,tol=tol)
    end

    return copy(boundaries)
end
