function leftorth(t::Array{T,N},leftind::NTuple{N1,Int},rightind::NTuple{N2,Int}) where {T,N,N1,N2}
    nt = permutedims(t,tuple(leftind...,rightind...))

    leftdims = broadcast(x->size(t,x),leftind)
    rightdims = broadcast(x->size(t,x),rightind)

    (Q,R) = qr(reshape(nt,prod(leftdims),prod(rightdims)))

    mid_dim = size(R,1);

    return (reshape(Matrix(Q),tuple(leftdims...,mid_dim)),reshape(Matrix(R),tuple(mid_dim,rightdims...)))
end

function rightorth(t::Array{T,N},leftind::NTuple{N1,Int},rightind::NTuple{N2,Int}) where {T,N,N1,N2}
    nt = permutedims(t,tuple(leftind...,rightind...))

    leftdims = broadcast(x->size(t,x),leftind)
    rightdims = broadcast(x->size(t,x),rightind)

    (L,Q) = lq(reshape(nt,prod(leftdims),prod(rightdims)))

    mid_dim = size(L,2);

    return (reshape(Matrix(L),tuple(leftdims...,mid_dim)),reshape(Matrix(Q),tuple(mid_dim,rightdims...)))
end

function Base.reverse(b :: Zygote.Buffer)

    bc = copy(b);
    tor = Zygote.Buffer(bc);

    for i = 1:length(bc)
        tor[i] = bc[end-i+1]
    end

    return tor
end

function init_envs(above,pepsline,below)
    T = eltype(above[1]);

    temp = eye(T,1,1)#Matrix{T}(I,1,1);

    ltemp = eye(T,size(pepsline[1],West),size(pepsline[1],West))#Matrix{T}(I,size(pepsline[1],West),size(pepsline[1],West))
    rtemp = eye(T,size(pepsline[length(pepsline)],East),size(pepsline[length(pepsline)],East))#Matrix{T}(I,size(pepsline[end],East),size(pepsline[end],East))

    @ein leftstart[-1,-2,-3,-4] := ltemp[-2,-3]*temp[-1,-4]
    @ein rightstart[-1,-2,-3,-4] := rtemp[-2,-3]*temp[-1,-4]

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

    temp = eye(T,1,1)#Matrix{T}(I,1,1);


    @ein leftstart[-1,-2] := temp[-1,-2]
    @ein rightstart[-1,-2] := temp[-1,-2]

    NL = Zygote.Buffer([leftstart]);NR = Zygote.Buffer([rightstart]);
    NL[1] = leftstart;NR[1] = rightstart;
    for i in 1:(length(below)-1)
        push!(NL,mps_apply_transfer_left(NL[i],below[i],below[i]));
        push!(NR,mps_apply_transfer_right(NR[i],below[length(below)-i+1],below[length(below)-i+1]));
    end
    NR = reverse(NR);

    return (NL,NR)
end

#qr free vomps
function north_vomps(above,pepsline,obelow;tol=Defaults.tol,maxiter=Defaults.maxiter)
    T = eltype(above[1]);

    below = Zygote.bufferfrom(obelow);

    (GL,GR) = init_envs(above,pepsline,below);
    (NL,NR) = init_envs(below);

    # accumulate totalnorm seperately, storing well normalized tensors in below
    # not sure if it's needed, but otherwise the first tensor has an absuurd norm, the others are normal
    # now I just evenly distribute them at the end of vomps
    totalnorm = 1.0;

    err = 0.0;
    for it = 1:maxiter

        err = 0.0;

        for i = 1:length(above)-1
            @ein temp[-1,-2,-3,-4]:= (GL[i])[-1,7,8,9]*(above[i])[9,5,3,1]*(GR[i])[1,4,2,-4]*(pepsline[i])[7,-2,4,5,6]*conj(pepsline[i])[8,-3,2,3,6]
            @ein temp[-1,-2,-3,-4]:= pinv(NL[i])[1,-1]*temp[1,-2,-3,2]*pinv(NR[i])[-4,2]

            @ein angle[]:=temp[1,2,3,4]*conj(below[i])[1,2,3,4]
            fid = (angle[1]*angle[1]')/(norm(temp)^2*norm(below[i])^2)
            err = max(err,1-abs(fid))

            totalnorm = norm(temp);
            below[i] = temp/totalnorm;

            GL[i+1] = mps_apply_transfer_left(GL[i],pepsline[i],above[i],below[i])
            NL[i+1] = mps_apply_transfer_left(NL[i],below[i],below[i])
        end

        for i = length(above):-1:2
            @ein temp[-1,-2,-3,-4]:= (GL[i])[-1,7,8,9]*(above[i])[9,5,3,1]*(GR[i])[1,4,2,-4]*(pepsline[i])[7,-2,4,5,6]*conj(pepsline[i])[8,-3,2,3,6]
            @ein temp[-1,-2,-3,-4]:= pinv(NL[i])[1,-1]*temp[1,-2,-3,2]*pinv(NR[i])[-4,2]

            @ein angle[]:=temp[1,2,3,4]*conj(below[i])[1,2,3,4]
            fid = (angle[1]*angle[1]')/(norm(temp)^2*norm(below[i])^2)
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

    #err > tol && @warn "vomps failed to converge $(err)"

    return copy(below)
end

#=
This is the actual vomps I typically use, which is inverse free (contains qr's)
=#
function _north_vomps(above,pepsline,obelow;tol=Defaults.tol,maxiter=Defaults.maxiter)
    T = eltype(above[1]);

    below = Zygote.bufferfrom(obelow);

    (GL,GR) = init_envs(above,pepsline,below);

    err = 0.0;
    for it = 1:maxiter

        err = 0.0;

        for i = 1:length(above)-1
            @ein temp[-1,-2,-3,-4]:= (GL[i])[-1,7,8,9]*(above[i])[9,5,3,1]*(GR[i])[1,4,2,-4]*(pepsline[i])[7,-2,4,5,6]*conj(pepsline[i])[8,-3,2,3,6]

            @ein angle[]:=temp[1,2,3,4]*conj(below[i])[1,2,3,4]
            fid = (angle[1]*angle[1]')/(norm(temp)^2*norm(below[i])^2)
            err = max(err,1-abs(fid))

            (below[i],c) = leftorth(temp,(1,2,3),(4,));

            @ein below[i+1][-1,-2,-3,-4]:=c[-1,1]*below[i+1][1,-2,-3,-4]

            GL[i+1] = mps_apply_transfer_left(GL[i],pepsline[i],above[i],below[i])
        end

        for i = length(above):-1:2
            @ein temp[-1,-2,-3,-4]:= (GL[i])[-1,7,8,9]*(above[i])[9,5,3,1]*(GR[i])[1,4,2,-4]*(pepsline[i])[7,-2,4,5,6]*conj(pepsline[i])[8,-3,2,3,6]

            @ein angle[]:=temp[1,2,3,4]*conj(below[i])[1,2,3,4]
            fid = (angle[1]*angle[1]')/(norm(temp)^2*norm(below[i])^2)
            err = max(err,1-abs(fid))


            (c,below[i]) = rightorth(temp,(1,),(2,3,4));

            @ein below[i-1][-1,-2,-3,-4]:=below[i-1][-1,-2,-3,1]*c[1,-4]

            GR[i-1] = mps_apply_transfer_right(GR[i],pepsline[i],above[i],below[i])
        end

        if err < tol
            break
        end

    end

    err > tol && @warn  "vomps failed to converge $(err)"

    return copy(below)
end
