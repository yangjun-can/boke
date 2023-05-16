%% 最大公因数

function [out]=undivide(a,b)

if mod(a,b)==0
    out=b;
else
    temp=undivide(b, mod(a,b));
    out=temp;
end

end