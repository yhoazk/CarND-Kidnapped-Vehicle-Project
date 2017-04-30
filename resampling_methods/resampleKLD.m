function [indx]=resampleKLD(q,S,KLD)
% [indx]=resampleKLD(q,S,KLD)
% KLD-resampling method for particle filtering; Current version is only for  
% the state as 1-dimension for simplicity. If multi-dimension is used, 
% replace S(j) with S(j,:) and correspondingly change the grid partitioning
% multi-dimension grid creating can learn from resampleDeterministic.m

% Author: Tiancheng Li,Ref:
% T. Li, S. Sun and T. Sattar, Adapting sample size in particle
% filters through KLD-resampling, Electronics Letters, vol. 49, no. 12, June
% 2013, p. 740-742
% based on KLD measure for resampling, the sample size is determined by
% the number of non-empty bins,k   

% Also see
% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering,
% submit to IEEE Signal Processing Magazine, August 2013

% Input:
%       q    the input weight sequence
%       S    the input state sequence 
%       KLD  KLD parameters£º
%            KLD.L is the bin size; KLD.Nmax is the largest sample size
% Output:
%       indx the resampled index according to the weight sequence
%% 
N = length(q);
xkld = ones(1,N) * 0.5;  
indx = zeros(1,N);
k = 0;  % initial number of non-empty bin
% nx=10; % we can set minimum sample size, if necessary 10 is because
% when k=2 (two non-empty bins, (k-1)/(2*KLD.e)*(1-2/(9*k-9)+sqrt(2/(9*k-9))*norminv(1-KLD.delta)^3)> 9 approximately )
%% Resampling start
% Based on Multinomial resamplng to select particles one by one
Q = cumsum(q);
Q(numel(q)) = 1;
i = 0;
while i < KLD.Nmax
  j = 1;
  simpl = rand;
  while Q(j) < simpl
    j = j + 1;
  end;
    i = i + 1;
    indx(i) = j;
    xkldpart = round((S(j)/KLD.L));
%KLD sampling
    flag = 0;
    % check whether new particles falls into new empty bin  
    % actually this can be runned by the commond of Matlab
    %   ismember(xkldpart,xkld);
    h = 0;
    while h < k
        h = h + 1;
        if xkldpart == xkld(h)
            flag = 1;
            break;                   
        end
    end
    if flag == 0
       % if new particle falls into empty bin, update: 
       % the non-empty bin set,number of non-empty bins, required number of particles nx 
       k = k + 1;
       xkld (k) = xkldpart;
       nx = (k-1)/(2*KLD.e)*(1-2/(9*k-9)+sqrt(2/(9*k-9))*norminv(1-KLD.delta))^3;               
    end 
    if i >= nx
        break; % the KL-bound is reached
    end 
end
indx(i+1:end) = [];
