using Revise,Zygote,BackwardsLinalg,OMEinsum,fpeps,LinearAlgebra

let
    width = 4;
    height = 4;

    D = 2;
    d = 2;
    chi = 2;

    peps = Array{Array{ComplexF64,5},2}(undef,width,height);
    for i = 1:width
        for j = 1:height
            peps[i,j] = rand(ComplexF64,D,D,D,D,d);
        end
    end

    sx = 0.5*ComplexF64[0 1;1 0]; sy = 0.5*ComplexF64[0 1im;-1im 0]; sz = 0.5*ComplexF64[1 0;0 -1];

    @ein ham1[-1,-2,-3,-4] := sx[-1,-2]*sx[-3,-4]
    @ein ham2[-1,-2,-3,-4] := sy[-1,-2]*sy[-3,-4]
    @ein ham3[-1,-2,-3,-4] := sz[-1,-2]*sz[-3,-4]

    ham = ham1+ham2+ham3;

    stepsize = 0.01;

    boundaries = gen_boundaries(peps,chi)
    for topit = 1:500


        (energy,boundaries) = calc_energy(peps,ham,boundaries)
        @show energy/(width*height);

        function cfun(x)
            calc_energy(x,ham,boundaries)[1]
        end

        derivs = cfun'(peps);
        @show sum(norm.(derivs))
        peps.-=stepsize*derivs;

    end

end
