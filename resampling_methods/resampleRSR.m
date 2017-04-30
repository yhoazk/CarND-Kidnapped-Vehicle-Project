function [ indx ] = resampleRSR(w, N)
% [ indx ] = resampleRSR(w, N)
% Residual Systematic resampling method for particle filtering.  
% Author: Tiancheng Li,Ref:
% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
% submit to IEEE Signal Processing Magazine, August 2013

% Input:
%       w    the input weight sequence 
%       N    the desired length of the output sequence(i.e. the desired number of resampled particles)
% Output:
%       indx the resampled index according to the weight sequence

if nargin == 1
  N = length(w);
end
M = length(w);
w = w / sum(w);
indx = zeros(1, N);

i =1;
u = rand/N;
j = 0;
while j < M
    j = j + 1;
    Ns = floor(N*(w(j)-u))+1;
      counter = 1;
      while counter <= Ns
        indx(i) = j;
        i = i + 1; counter = counter + 1;
      end;
    u = u + Ns/N-w(j);
end;

%% the following also works but may be  slow
% w = w ./ sum(w);
% i=0;
% u=rand/N;
% for j=1:N,
%     Ns=floor(N*(w(j)-u))+1;
%     for k=1:Ns
%         i = i +1;
%         indx(i)=j;
%     end
%     u=u+Ns/N-w(j);
% end;

%  The following lines are slow
%     if Ns
%         indx(i:i+Ns)=j;
%         i = i +Ns;
%     end
