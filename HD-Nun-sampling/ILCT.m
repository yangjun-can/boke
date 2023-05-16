function [x,ts] = ILCT(X_lct,a,b,c,d, M, us)
N=M;
S=sign(b);
% k1=sqrt(1/M);
k1=sqrt(1i*2*pi*b)/N;
u=us*(0:M-1);
chirp_u = exp(-1i*d/(2*b)*u.^2);
x_t=X_lct.*chirp_u;
if b>0
     x=ifft(x_t)*N; 
   else
      x=fft(x_t); 
end
ts=2*pi*S*b/(us*N);
t=ts*(0:N-1);
chirp_t = exp(-1i*a/(2*b)*t.^2);
x=k1*chirp_t.*x;
end