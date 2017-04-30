function [xpart,W,Num]=resampleDeterministic(N,w,x,DR)
% [xpart,W,Num]=resampleDeterministic(N,w,x,DR)
% [xpart,W,Num]=resampleDeterministic(N,w,x,DR)
% Deterministic resampling for particle filtering.

% Author: Tiancheng Li,Ref:
% T. Li, T. P. Sattar, and S. Sun, “Deterministic resampling: 
% Unbiased sampling to avoid sample impoverishment in particle filters,” 
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

w = w /sum(w);
a = length(w);
S = x + 100; % 100 is a number of big enough to elimate all negetive number

[Gt,d,Wt] = divi2(S,DR.Lstar,w,DR.density,DR.Lmin);
% N is specified
w = w*N;
qm = floor(w);
% % Use roughening to prevent sample impoverishment.
%        E = max(S')' - min(S')';
%        sigma = 0.2 * E * N^(-1); %-1/length(x)
temp = 0;
j = 0;
while j< a
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
for i = 1:k
    if d(i) ~= 0
       qqt = sum(Wt(j+1:j+d(i)));
       if qqt >= 0.1/N % this could be zero
           temp = temp+1;
           xpart(n+temp) = sum(Gt(j+1:j+d(i)).*Wt(j+1:j+d(i)))/qqt;
           W(n+temp) = qqt; 
       end
    end
    j = j + d(i);
end
Num = n + temp;
if temp > 0
    qsum = sum(W(n+1:Num));
else
    qsum = 0;
end
W(1:n) = (1-qsum)*ones(1,n)/n;

xpart = xpart - 100; 

xpart(Num+1:end) = [];
W(Num+1:end) = [];



function [Gt,Dt,Wt] = divi2(S,L,q,density,Lmin)

Nt = size(q,2);
i = 1;
while i <= Nt
    if q(i) < 0.02/Nt
        q(i) = [];
        S(i) = [];
        Nt = Nt - 1;
    else
      i = i + 1;  
    end
end

  imax = round(max((S)/L));
  imin = round(min((S)/L));
Gt = [];
Wt = [];
Dt = [];
for i = imin : imax
    d(i)  = 0;
end
p = 0;
while p < Nt
    p = p + 1;
    i = round(S(p)/L);
    d(i) = d(i) + 1;
    P(i,d(i)) = S(p);
    W(i,d(i)) = q(p);
end
i = imin;
while i < imax
    if d(i) > 0
        T = [];
        w = [];
        for p = 1 : d(i)
           T = [T,P(i,p)];
           w = [w,W(i,p)];
        end
        if d(i) >= density & L >= Lmin
           [Gthat,Dthat,Wthat] = divi2(T,L/2,w,density,Lmin);
           for p = 1 : d(i)
               P(i,p) = 0;
               W(i,p) = 0;
           end
           d(i) = 0;
           Dt = [Dt,Dthat];
           Gt = [Gt,Gthat];
           Wt = [Wt,Wthat];
        else
           Gt = [Gt,T];
           Wt = [Wt,w];
           Dt = [Dt,d(i)];
        end
    end
    i = i + 1;
end

    