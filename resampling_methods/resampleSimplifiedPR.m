function [ indx ] = resampleSimplifiedPR( w, T, N )
% [ indx ] = resampleSimplifiedPR( w, T, N )
% Simplified Partial resampling method for particle filtering.  
% Author: Tiancheng Li,Ref:
% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
% submit to IEEE Signal Processing Magazine, August 2013

% Input:
%       w    the input weight sequence 
%       N    the desired length of the output sequence(i.e. the desired number of resampled particles)
%       T    threshold to divide particles into two groups
% Output:
%       indx the resampled index according to the weight sequence

if nargin == 1
  N = length(w);
  T = 5*N;
end
if nargin == 2
  N = length(w);
end
M = length(w);
w = w / sum(w);
indx = zeros(1, N);

h = 0;
j = 0;
while j < M
    j = j + 1;
    if w(j) > 1/T
        h = h + 1;
        A(h) = j;
    end
end
r = 1;
i = 0;
while i < N
    i = i + 1;
    indx(i) = A(r);
    r = mod(r,h) + 1;
end
