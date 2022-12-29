function [A,Zm] = sxlp(s,p,w,varargin)
% A = sxlp(s,p,w)
% A   - filter coefficients
% s   - signal
% p   - prediction order
% w   - weighting parameter
% ... - additional parameters for weighting in case w is a handle to external function

if all(s==0)
  s = eps*randn(size(s));
end

if nargin<3
  w = p;
end

N = length(s);

if nargout>1
  Zm = zeros(p+1,N+p);
end

D = flipud(buffer([s;zeros(p,1)],p+1,p));

if ~isa(w,'function_handle')
  % AVS weighting
  ac = zeros(p+1,1);
  mem = (w-1)./w;
  Zp = zeros(p+1,1);
  for i1=1:N+p
    ac = mem.*ac + (1-mem).*(abs(D(1,i1))+abs(D(:,i1)));
    Z = max([ac [0;Zp(1:p)]],[],2);
    D(:,i1) = D(:,i1) .* Z;
    Zp = Z;
    if nargout>1
      Zm(:,i1) = Z;
    end
  end
else
  % external weighting function
  Zm = feval(w,s,varargin{:});
  for i1=2:N+p
    Zm(:,i1) = max([Zm(:,i1) [0;Zm(1:p,i1-1)]],[],2);
  end
  D = D.*Zm;
end

R = (D*D')/N;
A = R(2:end,2:end)\R(2:end,1);
A = [1;-A]';
