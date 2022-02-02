function A=FastJLmat_RFD(R, N, D)
% Fast JL matrix RFD where R is restriction, F is DFT, and D is a random
% diagonal matrix
% Refer to corollary 3.2 in S. Bamberger and F. Krahmer, "Optimal Fast Johnson-Lindenstrauss Embeddings for Large Data"

% A: m x N fast JL matrix

% R: m-tuple restriction vector containing the indices that should be kept.

% N: ambient dimension

% D: Rademacher sequence (1 x N) used to form the random diagonal

% Created by Ali Zare

m=length(R);
A=exp(-1i*2*pi/N*(R-1).'*(0:N-1)).*D/sqrt(m);
