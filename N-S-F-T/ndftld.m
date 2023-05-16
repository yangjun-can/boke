function [Xk]= ndftld(xn, N, fs)
% compute NnUnifon Discrete Fourier Transform 
% ------------------------------------------
% xn=N-Doint 1D inout sequence over 0<=n<=N-1(Uniform)
% Xk=1D NDFT coefficient arry (NonUniform)
% fs= nonuniform frequency vector(N-point)
% N= Length of DFT
% Usage: Xk=ndft1a(xn，N, fs)
n=0:1 : N-1; %Index for input data
Wn=exp(-1i*2*pi/N); % Twiddle factor
nk=fs .* n;%Creates an NxN matrix
DFTmtx=Wn .^nk;%DFT matrix (NXN)
Xk = (DFTmtx* xn.');%DFT coefficients (co1umm vector)
end