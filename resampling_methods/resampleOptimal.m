function [ W, indx ] = resampleOptimal(w, N)
% [indx ] = resampleStratified(w, N)
% Optimal resampling method for particle filtering.
%  
% Author: T. Li, 
% Ref: J.Carpenter, and P. Fearnhead,   
% On-line inference for hidden Markov models via particle filters
% Also ref: T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
% IEEE Signal Processing Magazine

% Input:
%       w    the input weight sequence 
%       N    the desired length of the output sequence(i.e. the desired number of resampled particles)
% Output:
%       W    the weights of resampled particles 
%       indx the resampled index according to the weight sequence

if nargin == 1
    N = length(w);
end
M = length(w);
w = w / sum(w);

if N < M
    % search for ct that satisfies sum_{n}(min(1, w_{n}/ct)) = N
    [w_s, w_ind] = sort(w);
    N_p=0;
    for ct = w_s
        N_p = N_p+1;
        if (M-N_p + sum(w_s(1:N_p))/ct)<= N
            break;
        end
    end
    ct = sum(w_s(1:N_p))/(N-M+N_p);
    
    % part 1 % Resample weights greter than ct; w_s are sorted
    indx = w_ind(N_p+1:M);
    W = w_s(N_p+1:M);
    % part 2 % Resample weights smaller than ct
    indx_part2 = resampleStratified( w_s(1:N_p), N-M+N_p);
    indx = [indx w_ind(indx_part2)];
    W = [W ones(1,N-M+N_p)*ct];
else 
    % otherwise just keep the same as original
    W = w;
    indx = 1:M;
end

