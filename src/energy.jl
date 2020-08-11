#need to figure out how to both return the energy and cache intermediaries
function calc_energy(peps,ham,prev_boundaries)

    boundaries = calc_boundaries(peps,prev_boundaries)
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

        temp = eye(T,1,1);

        NLtemp = eye(T,size(peps[r,1],West),size(peps[r,1],West))
        NRtemp = eye(T,size(peps[r,end],East),size(peps[r,end],East))

        @ein NL[-1,-2,-3,-4] := temp[-1,-4]*NLtemp[-2,-3]
        @ein NR[-1,-2,-3,-4] := temp[-1,-4]*NRtemp[-2,-3]

        NL = NL
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

        @ein h[] := GL[1,2,3,4]*NR[4,2,3,1]
        @ein n[] := NL[1,2,3,4]*NR[4,2,3,1]

        toret += sum(h)/sum(n)

    end

    return toret
end
