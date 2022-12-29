function [A,w] = wlp(s,p,M,k)
% A = wlp(s,p,m,k)
% weighted linear prediction
% A - filter coefficients
% s - signal
% p - prediction order
% m - STE window length (scalar) or weighting function (vector)
% k - STE window lag (default 1)

  % J.P. 120509, 180909, 100910
  
s = s + eps;

N = length(s);

if nargin<4
    k = 1;
end

if length(M)==1
    w = stew(s,M,p,k)+eps;
else
    w = M;
    if length(w)<N+p
        error('Weighting function must have at least N+p values.');
    end
end

wsr = sqrt(w(1:(N+p)));
Y = zeros(N+p,p+1);
for i1=0:p
    Y(:,i1+1) = [zeros(i1,1);s;zeros(p-i1,1)].*wsr;
end

R = (Y'*Y)/N;
A = inv(R(2:end,2:end))*R(2:end,1);
A = [1;-A]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = stew(x,m,p,k)
% STE weighting
w = conv([zeros(k,1);x(:)].^2,ones(m,1));
w = [w;zeros(p-m,1)];
