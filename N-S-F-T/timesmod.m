%% 重排前位置out,重排后位置x
function [out]=timesmod(x, a, n) 
% out=fix(rem(x*a,n));      % rem:取余函数；fix:向0靠近取整;
out=mod(x*a,n);    
end