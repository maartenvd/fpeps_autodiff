const Dir = UInt8;

const North =   0b100 :: Dir;
const East  =   0b011 :: Dir;
const South =   0b010 :: Dir;
const West  =   0b001 :: Dir;

const Dirs = (West,South,East,North);

TensorKit.space(t::AbstractTensorMap, d::Dir) = space(t,Int64(d))

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
        dnpeps = map(x->permute(x,(4,1,2,3),(5,)),npeps);
        return rotate_north(dnpeps,left(dir))
    end
end

function rotate_north(pepst::AbstractTensorMap,dir::Dir) where T
    if dir == North
        return pepst
    else

        pepst = permute(pepst,(4,1,2,3),(5,));
        return rotate_north(pepst,left(dir))
    end
end

function Base.reverse(b :: Zygote.Buffer)

    bc = copy(b);
    tor = Zygote.Buffer(bc);

    for i = 1:length(bc)
        tor[i] = bc[end-i+1]
    end

    return tor
end
