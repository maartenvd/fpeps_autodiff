#usual vomps
struct QR_vomps
    tol::Float64
    maxiter::Int
end

QR_vomps(;tol=Defaults.tol,maxiter=Defaults.maxiter) = QR_vomps(tol,maxiter);

#usual vomps, but we invert reduced density matrices
struct QR_free_vomps
    tol::Float64
    maxiter::Int
end

QR_free_vomps(;tol=Defaults.tol,maxiter=Defaults.maxiter) = QR_free_vomps(tol,maxiter);

#2site vomps
struct SVD_vomps
    tol::Float64
    trunc::Float64
    maxiter::Int
end
SVD_vomps(;tol=Defaults.tol,trunc=Defaults.trunc,maxiter=Defaults.maxiter) = SVD_vomps(tol,trunc,maxiter);

function approximate(above,pepsline,obelow,alg::QR_free_vomps)
    tol = alg.tol; maxiter = alg.maxiter;

    T = eltype(above[1]);

    below = Zygote.bufferfrom(obelow);

    (GL,GR) = init_envs(above,pepsline,below);
    (NL,NR) = init_envs(below);

    # accumulate totalnorm seperately, storing well normalized tensors in below
    # not sure if it's needed, but otherwise the first tensor has an absurd norm, the others are normal
    # now I just evenly distribute them at the end of vomps
    totalnorm = 1.0;

    err = 0.0;
    for it = 1:maxiter

        err = 0.0;

        for i = 1:length(above)-1
            @tensor temp[-1 -2 -3;-4]:= (GL[i])[-1,7,8,9]*(above[i])[9,5,3,1]*(GR[i])[1,4,2,-4]*(pepsline[i])[7,-2,4,5,6]*conj(pepsline[i][8,-3,2,3,6])
            @tensor temp[-1 -2 -3;-4]:= pinv(NL[i])[-1,1]*temp[1,-2,-3,2]*pinv(NR[i])[2,-4]

            angle = @tensor temp[1,2,3,4]*conj(below[i][1,2,3,4])
            fid = (angle*angle')/(norm(temp)^2*norm(below[i])^2)
            err = max(err,1-abs(fid))

            totalnorm = norm(temp);
            below[i] = temp/totalnorm;

            GL[i+1] = mps_apply_transfer_left(GL[i],pepsline[i],above[i],below[i])
            NL[i+1] = mps_apply_transfer_left(NL[i],below[i],below[i])
        end

        for i = length(above):-1:2
            @tensor temp[-1 -2 -3;-4]:= (GL[i])[-1,7,8,9]*(above[i])[9,5,3,1]*(GR[i])[1,4,2,-4]*(pepsline[i])[7,-2,4,5,6]*conj(pepsline[i][8,-3,2,3,6])
            @tensor temp[-1 -2 -3;-4]:= pinv(NL[i])[-1,1]*temp[1,-2,-3,2]*pinv(NR[i])[2,-4]

            angle = @tensor temp[1,2,3,4]*conj(below[i][1,2,3,4])
            fid = (angle*angle')/(norm(temp)^2*norm(below[i])^2)
            err = max(err,1-abs(fid))

            totalnorm = norm(temp);
            below[i] = temp/totalnorm;

            GR[i-1] = mps_apply_transfer_right(GR[i],pepsline[i],above[i],below[i])
            NR[i-1] = mps_apply_transfer_right(NR[i],below[i],below[i])
        end

        if err < tol
            break
        end
    end

    totalnorm = totalnorm^(1/length(below));
    for i in 1:length(below)
        below[i]*=totalnorm
    end

    #cannot autodiff through try-catch so @warn doesn't work -.-
    #err > tol && @warn "vomps failed to converge $(err)"

    return copy(below)
end

function approximate(above,pepsline,obelow,alg::QR_vomps)
    tol = alg.tol; maxiter = alg.maxiter;

    T = eltype(above[1]);

    below = Zygote.bufferfrom(obelow);

    (GL,GR) = init_envs(above,pepsline,below);

    err = 0.0;
    for it = 1:maxiter

        err = 0.0;

        for i = 1:length(above)-1
            @tensor temp[-1 -2 -3;-4]:= (GL[i])[-1,7,8,9]*(above[i])[9,5,3,1]*(GR[i])[1,4,2,-4]*(pepsline[i])[7,-2,4,5,6]*conj(pepsline[i][8,-3,2,3,6])
            err = max(err,norm(temp-below[i])/norm(below[i]))
            (below[i],c) = leftorth(temp,(1,2,3),(4,));

            @tensor below[i+1][-1 -2 -3;-4]:=c[-1,1]*below[i+1][1,-2,-3,-4]

            GL[i+1] = mps_apply_transfer_left(GL[i],pepsline[i],above[i],below[i])
        end

        for i = length(above):-1:2
            @tensor temp[-1 -2 -3;-4]:= (GL[i])[-1,7,8,9]*(above[i])[9,5,3,1]*(GR[i])[1,4,2,-4]*(pepsline[i])[7,-2,4,5,6]*conj(pepsline[i][8,-3,2,3,6])
            err = max(err,norm(temp-below[i])/norm(below[i]))

            (c,t) = rightorth(temp,(1,),(2,3,4));
            below[i] = permute(t,(1,2,3),(4,))
            @tensor below[i-1][-1 -2 -3;-4]:=below[i-1][-1,-2,-3,1]*c[1,-4]

            GR[i-1] = mps_apply_transfer_right(GR[i],pepsline[i],above[i],below[i])
        end


        if err < tol
            break
        end

    end

    #err > tol && @warn "vomps failed to converge $(err)"

    return copy(below)
end

function approximate(above,pepsline,obelow,alg::SVD_vomps)
    tol = alg.tol; maxiter = alg.maxiter; trunc = truncbelow(alg.trunc);

    T = eltype(above[1]);

    below = Zygote.bufferfrom(obelow);

    (GL,GR) = init_envs(above,pepsline,below);

    err = 0.0;
    for it = 1:maxiter

        err = 0.0;

        for i = 1:length(above)-2
            @tensor temp[-1,-2,-3,-4,-5,-6] := GL[i][-1,1,2,3]*above[i][3,4,5,6]*above[i+1][6,7,8,9]*GR[i+1][9,10,11,-6]*pepsline[i][1,-2,12,4,14]*conj(pepsline[i][2,-3,13,5,14])*pepsline[i+1][12,-4,10,7,15]*conj(pepsline[i+1][13,-5,11,8,15])
            (U,S,V) = tsvd(temp,(1,2,3),(4,5,6),trunc = trunc)

            @tensor derp1[-1 -2 -3;-4 -5 -6]:=below[i][-1,-2,-3,1]*below[i+1][1,-4,-5,-6]
            @tensor derp2[-1 -2 -3;-4 -5 -6]:=U[-1,-2,-3,1]*S[1,2]*V[2,-4,-5,-6]
            err = max(err,norm(derp1-derp2)/norm(temp));

            below[i] = U
            @tensor below[i+1][-1 -2 -3;-4]:=S[-1,1]*V[1,-2,-3,-4]

            GL[i+1] = mps_apply_transfer_left(GL[i],pepsline[i],above[i],below[i])
        end

        for i = length(above)-1:-1:1
            @tensor temp[-1,-2,-3,-4,-5,-6] := GL[i][-1,1,2,3]*above[i][3,4,5,6]*above[i+1][6,7,8,9]*GR[i+1][9,10,11,-6]*pepsline[i][1,-2,12,4,14]*conj(pepsline[i][2,-3,13,5,14])*pepsline[i+1][12,-4,10,7,15]*conj(pepsline[i+1][13,-5,11,8,15])
            (U,S,V) = tsvd(temp,(1,2,3),(4,5,6),trunc = trunc)

            @tensor derp1[-1 -2 -3;-4 -5 -6]:=below[i][-1,-2,-3,1]*below[i+1][1,-4,-5,-6]
            @tensor derp2[-1 -2 -3;-4 -5 -6]:=U[-1,-2,-3,1]*S[1,2]*V[2,-4,-5,-6]
            err = max(err,norm(derp1-derp2)/norm(temp));

            @tensor below[i][-1 -2 -3;-4] := U[-1,-2,-3,1]*S[1,-4]
            below[i+1] = permute(V,(1,2,3),(4,))

            GR[i] = mps_apply_transfer_right(GR[i+1],pepsline[i+1],above[i+1],below[i+1])
        end

        if err < tol
            break
        end

    end

    err > tol && @warn  "vomps failed to converge $(err)"

    return copy(below)
end


function init_envs(above,pepsline,below)
    T = eltype(above[1]);

    ou = oneunit(space(below[1],1));
    temp = isometry(Matrix{T},ou,ou);

    ltemp = isometry(Matrix{T},space(pepsline[1],West)',space(pepsline[1],West)');
    rtemp = isometry(Matrix{T},space(pepsline[length(pepsline)],East)',space(pepsline[length(pepsline)],East)')

    @tensor leftstart[-1 -2 -3;-4] := ltemp[-2,-3]*temp[-1,-4]
    @tensor rightstart[-1 -2 -3;-4] := rtemp[-2,-3]*temp[-1,-4]

    GL = Zygote.Buffer([leftstart]);GR = Zygote.Buffer([rightstart]);
    GL[1] = leftstart;GR[1] = rightstart;
    for i in 1:(length(above)-1)
        push!(GL,mps_apply_transfer_left(GL[i],pepsline[i],above[i],below[i]));
        push!(GR,mps_apply_transfer_right(GR[i],pepsline[length(pepsline)-i+1],above[length(pepsline)-i+1],below[length(pepsline)-i+1]));
    end
    GR = reverse(GR);

    return (GL,GR)
end

function init_envs(below)
    T = eltype(below[1]);

    ou = oneunit(space(below[1],1));
    temp = isometry(Matrix{T},ou,ou);


    @tensor leftstart[-1;-2] := temp[-1,-2]
    @tensor rightstart[-1;-2] := temp[-1,-2]

    NL = Zygote.Buffer([leftstart]);NR = Zygote.Buffer([rightstart]);
    NL[1] = leftstart;NR[1] = rightstart;
    for i in 1:(length(below)-1)
        push!(NL,mps_apply_transfer_left(NL[i],below[i],below[i]));
        push!(NR,mps_apply_transfer_right(NR[i],below[length(below)-i+1],below[length(below)-i+1]));
    end
    NR = reverse(NR);

    return (NL,NR)
end
