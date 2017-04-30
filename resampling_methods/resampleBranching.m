function [ indx ] = resampleBranching(w,N)
% [ indx ] = resampleBranching(w,N)
% Branching-kill method for resampling for particle filtering.  
% Author: Tiancheng Li,Ref:
% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
% submit to IEEE Signal Processing Magazine, August 2013

% Input:
%       w    the input weight sequence 
%       N    the desired length of the output sequence(i.e. the desired number of resampled particles)
% Output:
%       indx the resampled index according to the weight sequence

if nargin ==1
  N = length(w);
end
M = length(w);
w = w / sum(w);

Nw = N * w;
Nw_bulk = floor(Nw); 

p = Nw - Nw_bulk;
flag = rand(1, M) < p; % 0 -- kill, 1 -- branch
new_num = Nw_bulk + flag;
N2 = sum(new_num);
% if N2 > 1
    indx = zeros(1,N2);
% end

i = 1;
j = 0;
while j < M
  cnt = 1;j = j + 1;
  while cnt <= new_num(j)
    indx(i) = j;
    i = i + 1; cnt = cnt + 1;
  end;
end;

