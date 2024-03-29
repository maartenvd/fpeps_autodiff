#use an algorithm from optimkit to minimize the energy
function find_groundstate(state,ham,derivalg,retractalg,optimalg,pars)

    function objfun(x)
        (state,pars,energy) = x;

        cfun(t) = calc_energy(t,ham,pars,alg=derivalg)[1]
        derivs = cfun'(state);
        derivs ./= size(state,1)*size(state,2)
        return energy,derivs
    end

    function retract(x, cgr, α)
        (ostate,opars,oenergy) = x;

        @info "trying stepsize $α"
        flush(stdout)
        flush(stderr)

        nstate = ostate+α*cgr;

        (nenergy,npars) = calc_energy(nstate,ham,opars,alg=retractalg);
        (nenergy,npars) = calc_energy(nstate,ham,npars,alg=derivalg);

        return (nstate,npars,nenergy),cgr
    end

    function inner(x, v1, v2)
        tot = 0.0;

        for (p1,p2) in zip(v1,v2)
            tot += real(dot(v1,v2))
        end

        return tot
    end
    function transport!(v, xold, d, α, xnew)
        v
    end
    scale!(v, α) = v.*α
    add!(vdst, vsrc, α) = vdst+α.*vsrc

    #return optimtest(objfun, (state,pars,calc_energy(state,ham,pars)[1]), objfun((state,pars,calc_energy(state,ham,pars)[1]))[2]; alpha= -0.01:0.001:0.01,retract = retract, inner = inner)

    pars = calc_boundaries(state,pars,alg=retractalg)
    (en,pars) = calc_energy(state,ham,pars,alg=derivalg)

    (x,fx,gx,normgradhistory)=optimize(objfun,(state,pars,en),optimalg;
        retract = retract,
        inner = inner,
        transport! = transport!,
        scale! = scale!,
        add! = add!,
        isometrictransport = false)

    return (x[1],x[2],normgradhistory[end])

end
