using Revise,TensorKit,TensorKitAD,fpeps,OptimKit

let
    width = 10;
    height = 10;

    D = ℂ^2;
    d = ℂ^2;
    chi = ℂ^10;

    peps = map(Iterators.product(1:width,1:height)) do (i,j)
        TensorMap(randn,ComplexF64,D*D*(D)'*(D)',(d)')
    end;

    sx = 0.5*TensorMap(ComplexF64[0 1;1 0],ℂ^2,ℂ^2);
    sy = 0.5*TensorMap(ComplexF64[0 1im;-1im 0],ℂ^2,ℂ^2);
    sz = 0.5*TensorMap(ComplexF64[1 0;0 -1],ℂ^2,ℂ^2);

    @tensor ham1[-1,-2,-3,-4] := sx[-2;-1]*sx[-4;-3];
    @tensor ham2[-1,-2,-3,-4] := sy[-2;-1]*sy[-4;-3];
    @tensor ham3[-1,-2,-3,-4] := sz[-2;-1]*sz[-4;-3];

    ham = ham1+ham2+ham3;

    retval = find_groundstate(peps,ham,QR_vomps(maxiter=10),QR_vomps(maxiter=10),ConjugateGradient(verbosity=2),gen_boundaries(peps,chi));
end
