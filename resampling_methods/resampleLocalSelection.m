function [ W, xpart] = resampleLocalSelection( w, x)
% [ W, xpart] = resampleLocalSelection( w, x)
% Local Selection (parallel) resampling method for particle filtering, 
% Note the state form: x(:,j), j is the index
% This method is designed for parallel processing but here IT runs in serial. Author: T. Li, Ref:
% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
% submit to IEEE Signal Processing Magazine, August 2013

% Input:
%       w    the input weight sequence 
% Output:
%       indx the resampled index according to the weight sequence

if nargin <= 1
  disp('too few inputs (both weight and state are required)');
end
M = length(w);
W = zeros(1, M);
xpart = x;

    T = w(M) + w(1) + w(2);
    T1 = w(M) / T;
    T2 = (w(M) + w(1)) / T;
    u = rand;
    if u <= T1
        W(1) = w(M);
        xpart(:,1) = x(:,M);
    elseif u <= T2
        W(1) = w(1);
        xpart(:,1) = x(:,1);
    else
        W(1) = w(2);
        xpart(:,1) = x(:,2);
    end
j = 1;
while j < M-1
    j = j + 1;
    T = w(j-1) + w(j) + w(j+1);
    T1 = w(j-1) / T;
    T2 = (w(j-1) + w(j)) / T;
    u = rand;
    if u <= T1
        W(j) = w(j - 1);
        xpart(:,j) = x(:,j - 1);
    elseif u <= T2
        W(j) = w(j);
        xpart(:,j) = x(:,j);
    else
        W(j) = w(j + 1);
        xpart(:,j) = x(:,j + 1);
    end
end

    T = w(M-1) + w(M) + w(1);
    T1 = w(M-1) / T;
    T2 = (w(M-1) + w(M)) / T;
    u = rand;
    if u <= T1
        W(M) = w(M - 1);
        xpart(:,M) = x(:,M - 1);
    elseif u <= T2
        W(M) = w(M);
        xpart(:,M) = x(:,M);
    else
        W(M) = w(1);
        xpart(:,M) = x(:,1);
    end
