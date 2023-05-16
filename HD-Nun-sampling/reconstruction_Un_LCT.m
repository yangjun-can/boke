function ReX_Udis_lct=reconstruction_Un_LCT(X_NUdis_lct,a,b,c,d,us,ts,kesai,Non_uniform_sample_time,N)
InvA=exp(-1i*d/(2*b)*us^2*(0:N-1).^2);
InvC=exp(1i*d/(2*b)*us^2*(0:N-1).^2);
B1=repmat(exp(-1i*a/(2*b)*(2*(0:N-1).*kesai+kesai.^2)*ts^2).',1,N);
B2=exp(1i*2*pi/N*Non_uniform_sample_time.'*(0:N-1)); 
B=B1.*B2;
InvB = pinv(fft(B));
ReX_Udis_lct=N*InvC.*(InvB*(InvA.*X_NUdis_lct).').';
end