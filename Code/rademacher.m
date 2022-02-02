function A=rademacher(varargin)

% Rademacher random array
% Created by Ali Zare

N=nargin;

I=zeros(1,N);
for k=1:N
    I(k)=varargin{k};
end

A=2*(rand(I)<0.5) - 1;