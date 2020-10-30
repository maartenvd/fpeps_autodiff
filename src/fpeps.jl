module fpeps
    using Zygote,BackwardsLinalg,OMEinsum,LinearAlgebra,OptimKit

    #default settings
    module Defaults
        const maxiter = 10
        const tol = 1e-12
    end

    LinearAlgebra.default_svd_alg(A::Matrix) = LinearAlgebra.QRIteration();
    
    export calc_energy,gen_boundaries,calc_boundaries,find_groundstate

    #zygote helper?
    eye(T,p1,p2) = Matrix{T}(I,p1,p2)
    Zygote.@nograd eye

    include("rotations.jl")
    include("transfers.jl")
    include("vomps.jl")
    include("boundaries.jl")
    include("energy.jl")
    include("optimhook.jl")
end
