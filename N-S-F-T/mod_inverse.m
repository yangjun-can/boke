%% 求模逆 a*out mod n = 1
function [out]=mod_inverse(a, n)
i = n;
out = 0;
d = 1;
while a>0
    t = fix(i/a);
    x = a;
    a = mod(i,x);
    i = x;
    x = d;
    d = out - t*x;
    out = x;
end
out=rem(out,n);         % rem(x,y):求整除x/y的余数
if (out<0)
    out = rem((out+n),n);  
end
end