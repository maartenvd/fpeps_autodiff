using OMEinsum,fpeps,LinearAlgebra,OptimKit

let
    width = 6;
    height = 6;

    D = 2;
    d = 2;
    chi = 10;

    peps = Array{Array{ComplexF64,5},2}(undef,width,height);
    for i = 1:width, j = 1:height
        peps[i,j] = rand(ComplexF64,D,D,D,D,d);
        peps[i,j] /= norm(peps[i,j])
    end

    sx = 0.5*ComplexF64[0 1;1 0]; sy = 0.5*ComplexF64[0 1im;-1im 0]; sz = 0.5*ComplexF64[1 0;0 -1];

    @ein ham1[-1,-2,-3,-4] := sx[-1,-2]*sx[-3,-4]
    @ein ham2[-1,-2,-3,-4] := sy[-1,-2]*sy[-3,-4]
    @ein ham3[-1,-2,-3,-4] := sz[-1,-2]*sz[-3,-4]

    ham = ham1+ham2+ham3;

    retval = find_groundstate(peps,ham,ConjugateGradient(verbosity=2),gen_boundaries(peps,chi))
end
