#generate inital north boundary mpses of bond dimension chi
function gen_north_bounds(peps,chi)
    T = eltype(peps[1,1])

    boundaries = map(1:size(peps,1)) do i
        lds = [oneunit(chi)]
        rds = [oneunit(chi)]

        for j = 1:size(peps,2)
            push!(lds,infimum(fuse(space(peps[i,j],North)'*space(peps[i,j],North)*lds[j]),chi))
            push!(rds,infimum(fuse(space(peps[i,end-j+1],North)'*space(peps[i,end-j+1],North)*rds[j]),chi))
        end
        rds = reverse(rds)

        cb = map(1:size(peps,2)) do j
            psize = space(peps[i,j],North)

            if i>1
                lD = infimum(lds[j],rds[j])
                rD = infimum(lds[j+1],rds[j+1])
                TensorMap(randn,T,lD*psize'*psize,rD)
            else
                lD = oneunit(chi)
                rD = oneunit(chi)
                permute(isometry(Matrix{T},lD*psize',rD*psize'),(1,2,4),(3,))
            end

        end
    end
end

function gen_boundaries(peps,chi = 1)
    T = eltype(peps[1,1]);

    bs = map([Dirs...]) do dir
        gen_north_bounds(rotate_north(peps,dir),chi)
    end
end

function calc_north_boundaries(peps,pboundaries;alg=QR_free_vomps())
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

        boundaries[i+1] = approximate(boundaries[i],pepsline,boundaries[i+1],alg)

    end
    return copy(boundaries)
end

function calc_boundaries(peps,old_boundaries;alg=QR_free_vomps())
    boundaries = Zygote.Buffer(old_boundaries)

    for cd in Dirs
        boundaries[cd] = calc_north_boundaries(rotate_north(peps,cd),old_boundaries[cd];alg=alg)
    end

    return copy(boundaries)
end
