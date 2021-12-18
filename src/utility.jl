const Dir = UInt8;

const North =   0b100 :: Dir;
const East  =   0b011 :: Dir;
const South =   0b010 :: Dir;
const West  =   0b001 :: Dir;

const Dirs = (West,South,East,North);

left(dir::Dir) = Dir(mod1(dir+1,4))
right(dir::Dir) = Dir(mod1(dir-1,4))

function brotl90(b)
    c = Zygote.Buffer(b);

    for i = 1:size(b,1)
        for j = 1:size(b,2)
            c[size(b,1)-j+1,i] = b[i,j];
        end
    end

    return copy(c);
end

function rotate_north(peps::Array{T,2},dir::Dir) where T
    if dir == North
        return peps
    else
        npeps = brotl90(peps);
        dnpeps = map(x->permutedims(x,[2,3,4,1,5]),npeps);
        return rotate_north(dnpeps,left(dir))
    end
end

function rotate_north(pepst::Array{T,5},dir::Dir) where T
    if dir == North
        return pepst
    else

        pepst = permutedims(pepst,[2,3,4,1,5]);
        return rotate_north(pepst,left(dir))
    end
end


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

function tsvd(t::Array{T,N},leftind::NTuple{N1,Int},rightind::NTuple{N2,Int};trunc=0) where {T,N,N1,N2}
    nt = permutedims(t,tuple(leftind...,rightind...))

    leftdims = broadcast(x->size(t,x),leftind)
    rightdims = broadcast(x->size(t,x),rightind)

    (U,S,Vt) = svd(reshape(nt,prod(leftdims),prod(rightdims)));

    fi = findfirst(x->x<trunc*S[1],S);
    if isnothing(fi)
        mid_dim = length(S);
    else
        mid_dim = fi;
    end

    tU = U[:,1:mid_dim];
    tS = diagm(S[1:mid_dim]);
    tV = Vt'[1:mid_dim,:]

    return (reshape(tU,tuple(leftdims...,mid_dim)),tS,reshape(Matrix(tV),tuple(mid_dim,rightdims...)))
end
function Base.reverse(b :: Zygote.Buffer)

    bc = copy(b);
    tor = Zygote.Buffer(bc);

    for i = 1:length(bc)
        tor[i] = bc[end-i+1]
    end

    return tor
end
