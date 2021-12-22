#need to figure out how to both return the energy and cache intermediaries
function calc_energy(peps,ham,prev_boundaries;alg=QR_vomps())

    boundaries = calc_boundaries(peps,prev_boundaries,alg=alg)
    toret = 0.0 + 0im;

    toret += calc_energy_impl(peps,boundaries[North],boundaries[South],ham)
    toret += calc_energy_impl(rotate_north(peps,East),boundaries[East],boundaries[West],ham)

    return (real(toret),boundaries)
end


function calc_energy_impl(peps,nbounds,sbounds,ham)
    T = eltype(peps[1,1]);

    toret = 0.0 + 0im;

    for r = 1:size(peps,1)
        n = nbounds[r];
        s = sbounds[size(peps,1)-r+1];

        ou = oneunit(space(n[1],1));
        temp = isometry(Matrix{T},ou,ou);

        NLtemp = isometry(Matrix{T},space(peps[r,1],West),space(peps[r,1],West))
        NRtemp = isometry(Matrix{T},space(peps[r,end],East),space(peps[r,end],East))

        @tensor NL[-1,-2,-3,-4] := temp[-1,-4]*NLtemp[-3,-2]
        @tensor NR[-1,-2,-3,-4] := temp[-1,-4]*NRtemp[-3,-2]

        GL = crosstransfer(NL,peps[r,1],n[1],s[end],dir=West)*0;

        for c in 1:size(peps,2)-1
            GL = crosstransfer(GL,peps[r,c+1],n[c+1],s[end-c],dir=West)

            tt = hamtransfer(s[end-c],s[end-c+1],n[c],n[c+1],NL,
                rotate_north(peps[r,c],West),rotate_north(peps[r,c+1],West),ham)

            GL += hamtransfer(s[end-c],s[end-c+1],n[c],n[c+1],NL,
                rotate_north(peps[r,c],West),rotate_north(peps[r,c+1],West),ham)
            NL = crosstransfer(NL,peps[r,c],n[c],s[end-c+1],dir=West)
        end

        NL = crosstransfer(NL,peps[r,end],n[end],s[1],dir=West)

        h = @tensor GL[1,2,3,4]*NR[4,2,3,1]
        n = @tensor NL[1,2,3,4]*NR[4,2,3,1]
        
        toret += h/n

    end

    return toret
end
