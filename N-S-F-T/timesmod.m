%% ����ǰλ��out,���ź�λ��x
function [out]=timesmod(x, a, n) 
% out=fix(rem(x*a,n));      % rem:ȡ�ຯ����fix:��0����ȡ��;
out=mod(x*a,n);    
end