function crosstransfer(v,o,a,b;dir=North,bo=o)
    o = rotate_north(o,dir);
    bo = rotate_north(bo,dir);

    @ein v[-1,-2,-3,-4]:=v[1,2,3,4]*a[4,5,6,-4]*o[7,-2,5,2,8]*conj(bo)[9,-3,6,3,8]*b[-1,7,9,1]
end

#used for fixpoint equations
function mps_apply_transfer_left(v,o,a, b,bo=o)
    @ein v[-1,-2,-3,-4]:=v[1,2,3,4]*a[4,5,6,-4]*o[2,7,-2,5,8]*conj(bo)[3,9,-3,6,8]*conj(b)[1,7,9,-1]
end

function mps_apply_transfer_right(v,o,a, b,bo=o)
    @ein v[-1,-2,-3,-4]:=v[1,2,3,4]*a[-1,5,6,1]*o[-2,7,2,5,8]*conj(bo)[-3,9,3,6,8]*conj(b)[-4,7,9,4]
end

#nn thing
#m1;m2 = left
#cbt = top tensor
#m3;m4 = right
#[t1;t2]
function hamtransfer(m1,m2,m3,m4,cbt,t1,t2,nn;bt1=t1,bt2=t2)
@ein toret[-1,-2,-3,-4]:=m1[-1,18,19,20]*m2[20,4,2,1]*cbt[1,5,3,6]*m3[6,7,8,11]*m4[11,12,15,-4]*
t1[4,13,7,5,9]*t2[18,-2,12,13,14]*conj(bt1)[2,16,8,3,10]*conj(bt2)[19,-3,15,16,17]*
nn[9,10,14,17]
return toret
end
