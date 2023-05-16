%% ÇÐ±ÈÑ©·ò´°º¯Êý
function [out]= Cheb ( m, x)
pd=abs(x) <= 1;
out(pd)=cos(m * acos(x(pd)));
out(~pd)=real(cosh(m * acosh(x(~pd))));
% if abs(x) <= 1
%     out=cos(m * acos(x)); % acos(x):·´ÓàÏÒarccos(x) 
% else
%     out=real(cosh(m * acosh(x))); % cosh(x):Ë«ÇúÓàÏÒº¯Êý=(e^x+e^-x)/2
% end
end