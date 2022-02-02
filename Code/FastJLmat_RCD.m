function A=FastJLmat_RCD(R, N, D)
% Fast JL with DCT
% Fast JL matrix RFD where R is restriction, F is DCT type 1, and D is a random
% diagonal matrix
% Check FastJLmat_RCD.m for comparison.

% A: m x N fast JL matrix

% R: m-tuple restriction vector containing the indices that should be kept.

% N: ambient dimension

% D: Rademacher sequence (1 x N) used to form the random diagonal

% Created by Ali Zare

H=DCT_1(N);
H=H(R,:);
m=length(R);
A=H*diag(D)/sqrt(m);
