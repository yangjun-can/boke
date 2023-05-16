function deltaN=Periodic_sinc(para,K,n)
if para/K==fix(para/K)
    deltaN=1;
else
    deltaN=sin(pi*para*n/K)/n/sin(pi*para/K);
end