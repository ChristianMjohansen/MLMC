% https://people.maths.ox.ac.uk/gilesm/
% calculation of Brownian Bridge
%
% function dW = bb(Z,T)
%
% T       time interval
% Z(:,:)  input vectors of unit Normal variables

function W = bb(W,T)
  K = 1;
  N = size(W,1);

%
% perform Brownian Bridge interpolation
%

M = round(log(N)/log(2));
if N~=2^M
  error('error: number of timesteps not a power of 2')
end

for m = 1:M
  W([1:2:2^m 2:2:2^m],:) = ...
   [ 0.5*W(1:2^(m-1),:)+W(2^(m-1)+1:2^m,:)/sqrt(2)^(m+1) ; ...
     0.5*W(1:2^(m-1),:)-W(2^(m-1)+1:2^m,:)/sqrt(2)^(m+1) ];
end

W = sqrt(T)*W;

return