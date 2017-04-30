function [ indx ] = resampleRounding( w, N )
% [ indx ] = resampleRounding( w, N )
% Rounding-copy resampling for particle filtering.  
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
w = w / sum(w);
M = length(w);
Ns = round(N * w);
N2 = sum(Ns);
% if N2 > 1
    indx = zeros(1,N2);
% end
i = 0;
j = 0;
while j < M
  cnt = 1; j= j+1;
  while cnt <= Ns(j)
    i = i + 1; cnt = cnt + 1;  
    indx(i) = j;
  end; 
end;

%% optional choice 
% the following short codes are too slow
% %         indx(i:i+Ns(j))=j;
% %         i = i +Ns(j);

% % the following is also not as good as the one used above
% % --- computing speed trick in Matlab
% % it seems for ... end is computing slower than while... end

% if nargin == 1
%   N = length(w);
% end
% w = w / sum(w);
% M = length(w);
% Ns = round(N * w);
% N2 = sum(Ns);
% % if N2 > 1
%     indx = zeros(1,N2);
% % end
% i = 0;
% for j = 1 : M
%   counter = 1;
%   while counter <= Ns(j)
%     i = i + 1; counter = counter + 1;  
%     indx(i) = j;
%   end;    
% end;
