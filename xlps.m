function A = xlps(s,p,weighting,smoothing)
% A = XLPS(S,P,WEIGHTING,SMOOTHING)
% Extended weighted linear prediction using the autocorrelation snapshot
% (Pohjalainen,Alku 2013)
% Outputs:
% A - all-pole filter coefficients
% Inputs
% S - a signal frame
% P - prediction order
% WEIGHTING - a code indicating the weighting scheme to use, 1 and 2
% correspond to the article
% SMOOTHING - nonzero to use spectral smoothing operation (Section 2.3), 0
% for no smoothing

if nargin<4
    if weighting==1
        smoothing = true;
    else
        smoothing = false;
    end
end

s = s(:,1);

N = length(s);

b = (p-1)/p;
idx = (p+1):-1:1;
Q = zeros(p+1,p+1);
Qprev = Q;
R = zeros(p+1,p+1);

s = [zeros(p,1);s;zeros(p,1)];

X = s*s';
Y_1 = abs(s)*ones(1,N+2*p);
Y_1 = Y_1+Y_1';

switch weighting
    case 1
        Q = abs(s(p+1)) + Y_1(idx,idx);
    case 2
        Q = s(p+1).^2 + abs(X(idx,idx));
end

for i1=1:N+p
    switch weighting
        case 1
            Q = b*Q + (1-b)*(abs(s(i1+p)) + Y_1(i1+idx-1,i1+idx-1));
        case 2
            Q = b*Q + (1-b)*(s(i1+p).^2 + abs(X(i1+idx-1,i1+idx-1)));            
    end
    if smoothing
        Q = max(Qprev,Q);
        Qprev = [zeros(1,p+1);[zeros(p,1) Q(1:end-1,1:end-1)]];
    end
    R = R + X(i1+idx-1,i1+idx-1).*Q;
end

A = R(2:end,2:end)\R(2:end,1);
A = [1;-A]';