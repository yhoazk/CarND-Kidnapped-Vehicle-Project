function [xpart,W,Num] = resampleDeterministic1F(N,w,x,DR)
% [xpart,W,Num] = resampleDeterministic1F(N,w,x,DR)
% Deterministic resampling using Fixed grid size in one dimension
% 1F: 1-deminsional Fixed grid size

% Author: Tiancheng Li,Ref:
% T. Li, T. P. Sattar, and S. Sun, “Deterministic resampling: 
% Unbiased sampling to avoid sample impoverishment in particle filters,?
% Signal Processing, vol. 92(7), pp. 1637-1645, 2012.

% Also see reference
% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering,
% submit to IEEE Signal Processing Magazine, August 2013

% Input:
%       N    the reference number of particles,i.e. desired sample size
%       w    the input weight sequence
%       x    the input state sequence 
%       DR   Dr parameters for adaptive grid division
%            including three elements: L,density,Lmin
% Output:
%       W    the output weight sequence
%       xpart the output state sequence 
%       Num   the output number of particles
if nargin < 4
  N = length(w);
end
% allocate memory:
xpart = zeros(1,N);
W = zeros(1,N);

w = w/sum(w);
M = length(w);
S = x + abs(min(x)) + 1; % to deal with negative numbers

[Gt,d,Wt] = Divid(S,DR,w); % need state space division to creat grids
% Gt: particles stored in grids
% d: the number of particles in each grid
% Wt: weight stored in grids in order consistent with Gt

w = w*N;
qm = floor(w);
% 1stage: residual resampling 
temp = 0;
j = 0;
while j < M
    j = j + 1;
    cnt = 1;
    i = temp;
    while cnt <= qm(j)
       i = i + 1;cnt = cnt + 1;
       xpart(i) = S(j);
    end;
       temp = temp + qm(j);
       w(j) = w(j) - qm(j);
end
n = sum(qm);
Wt = (Wt*N-floor(Wt*N))/N;

temp = 0;
j = 0;
k = size(d,2);
i = 0;
while i < k
    i = i + 1;
    if d(i)~= 0
       qqt = sum(Wt(j+1:j+d(i)));
       if qqt > 0  % this could be non-zero so that to discard the over-small weighted particles  such as 0.1/N
           temp = temp+1;
           xpart(n+temp) = sum(Gt(j+1:j+d(i)).*Wt(j+1:j+d(i)))/qqt;
           W(n+temp) = qqt; 
       end
    end
    j = j + d(i);
end
Num = n + temp;
if temp>0
    qsum  =  sum(W(n+1:Num));
else
    qsum = 0;
end
W(1:n) = (1-qsum)*ones(1,n)/n;
xpart = xpart-abs(min(x))-1; 
xpart(Num+1:end) = [];
W(Num+1:end) = [];



function [Gt,Dt,Wt] = Divid(S,L,q)

Nt = size(S,2);
  imax = round(max((S)/L));
  imin = round(min((S)/L));
Gt = [];
Wt = [];
Dt = [];
    d(imin:imax)  = 0;
p = 0;
while p < Nt
    p = p + 1;
    i  =  round(S(p)/L);
    d(i) =  d(i) + 1;
    P(i,d(i)) = S(p);
    W(i,d(i)) = q(p);
end
i = imin;
while i < imax
    if d(i) > 0
        p = 0;
        while p < d(i)
            p = p + 1;
            Gt = [Gt,P(i,p)];
            Wt = [Wt,W(i,p)];
        end
           Dt = [Dt,d(i)];
    end
    i = i + 1;
end

    