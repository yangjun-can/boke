function [x,HSx,WSx]=signal(t)
% t : samples
% x : discreate signal
w0=0.005*pi;
w2=w0;
w1=0.1;
x = exp(-1i*w0*t.^2-w1*t) .*cos(w2*t.^2+w2*t) ;%.*heaviside(t)
%% Half-sample symmetric
sy_x=fliplr(x);%左右对称
HSx=cat(2,sy_x,x);  % 按列将x,syx连接起来
%% whole-sample symmetric
WSx=HSx;
WSx(length(t)+1)=[];
end