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
