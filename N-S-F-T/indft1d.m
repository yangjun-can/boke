function [xn] = indft1d(xk, N, fs)
%  Computes 1D NonUniform Inverse Discrete Fourier Transform
% 8-----------------------------------------
% xn = N-point 1D input sequence over 0 <= n <= N-1
% fs = Nonuniform frequency vector (N-point)
% Xk= 1D DFT coefficient array over 0 <= k<= N-1
% N= Length of DFT
% Usage: xn= indft1a(Xk, N, fs)
n=0: 1: N-1; % Index for input data
Wn= exp(-1i*2*pi/N); % Twiddle factor
nk= fs * n; % Createsan N X Nmatrix
DFTmtx= Wn.^nk; % DFT matrix (NX N)
xn= (inv(DFTmtx) *xk.') .';% Reconstructed sequence (row vector)
end