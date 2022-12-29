function [A,xi,B,C,Ps] = gmlp(s,p,Bi,iterations,postproblimit,lpcond,variancelimit)
% GMLP Gaussian mixture linear prediction
%        [A,xi,B,C,Ps] = gmlp(s,p,Bi,iterations,postproblimit,lpcond,variancelimit)
%
% OUTPUTS:
% A   - linear prediction inverse filter solved according to GMLP
% xi  - final posterior probabilities of different states after EM iteration
% Other outputs return parameter estimates after final EM iteration:
% B   - contains the autoregression parameters for each state in
% format explained below
% C   - variance parameters
% Ps  - component weights (prior probabilities)
%
% INPUTS:
% s  - signal frame
% p  - prediction order
% Bi - initial parameters in format [state1intercept state1predictorcoeffients;
% ...;stateJintercept stateJpredictorcoefficients] or "1" for GMLP-0 and "2" for GMLP-H
%      (DEFAULT 2)
% iterations - number of EM iterations (DEFAULT 3), if this is a vector then
% iterations(2) gives the number of preliminary iterations during which
% only the probability parameters (Ps) and variances (C) are updated
% (DEFAULT 0)
% postproblimit - a lower limit for posterior probabilities
%       (DEFAULT 0.01)
% lpcond - determines whether AR coefficients are solved according to the
%          autocorrelation method (solution conditional to lpcond=0 values)
%          or covariance method (solution conditional to lpcond=p values)
%       (DEFAULT 0)
% variancelimit - a lower limit for the variances
%       (DEFAULT 0.0)
%

% Jouni Pohjalainen Aalto University 18.11.2013 / 22.5.2014

if nargin<3
    Bi = 2;
end

if nargin<4
    iterations = 3;
end

if nargin<5
    postproblimit = 0.01;
end

if nargin<6
    lpcond = 0;
end

if nargin<7
    variancelimit = 0.0;
end

% boost frames that have a very low maximum absolute value
% in order to avoid underflow in computing autocorrelations
abslim = 10*sqrt(eps);
if all(abs(s)<abslim)
    s = s + abslim*randn(size(s));
end

N = length(s);

if isscalar(Bi)
    switch Bi
        case 1
            B = [0 0.97 zeros(1,p-1);zeros(1,p+1)];
        case 2
            B = [0 0.97 zeros(1,p-1);1 zeros(1,p)];
            iterations = [iterations 1];
    end
else
    B = Bi;
end

intercept = size(B,2)-p;

if intercept
    B(:,1:intercept) = B(:,1:intercept) * max(s);
end

xi = [];

if isscalar(iterations)
    iterations = [iterations 0];
end

J = size(B,1);
Ps = ones(J,1)/J;

C = repmat(0.01,1,J);

D = [flipud(buffer([0;s;zeros(p,1)],p,p-1))];
if intercept
    D = [ones(intercept,N+p+1);D];
end
N2 = N+p-lpcond;
s = [s;zeros(p-lpcond,1)];

for icount=0:(iterations(1)+iterations(2)-1)

    % E step

    res = repmat(s((lpcond+1):(N2+0))',J,1) - (B*D(:,(lpcond+1):N2)); % prediction residual
    Nr = size(res,2);
    exp_res = exp(-0.5.*repmat((1./C)',1,N2-lpcond).*(res.^2));
    % f(s(n)|j,s(n-1),...,s(n-p)), j=1,...,J
    eta = (1./sqrt((2*pi).*repmat(C',1,N2-lpcond))) .* exp_res + 1e-100;

    % f(s(n),j|s(n-1),...,s(n-p))
    xi = repmat(Ps,1,N2-lpcond) .* eta;
    % P(j|s(n),s(n-1),...,s(n-p))
    xi = xi ./ repmat(sum(xi,1),J,1);

    % enforce lower limit for the posterior probabilities
    xi = min(1-postproblimit,max(postproblimit,xi));
    % renormalize so that they sum to 1 for each observation
    xi = xi ./ repmat(sum(xi,1),J,1);

    % M step
    
    % update probability parameters
    Ps = sum(xi,2)./(N2-lpcond);

    C = zeros(size(C));
    cnorm = zeros(size(C));
    for j=1:J

        if icount>=iterations(2)
            % re-estimate autoregression parameters
            W  = repmat(sqrt(xi(j,:)),p+intercept,1);
            Zw = D(:,(lpcond+1):N2).*W;
            Xw = s((lpcond+1):(N2+0)).*sqrt(xi(j,:)');
            B(j,:) = (Xw'*Zw') / (Zw*Zw');

            % this would use the just computed new AR coefficients for residual
            % in order to re-estimate variance, similarly to Reynolds 1995
            % where updated mean values are used for variance
            %wres = (Xw - (B(j,:)*Zw)');
        end
        % this uses the previous iteration's AR values and the residual
        % already computed during the E step
        wres = res(j,:) .* sqrt(xi(j,:));
        whichcm = min([j length(C)]);
        C(whichcm) = C(whichcm) + sum(wres.^2,2);
        cnorm(whichcm) = cnorm(whichcm) + sum(xi(j,:));
    end

    % recompute variance and apply lower limit
    if length(variancelimit)>1
        C = min(variancelimit(2),max(variancelimit(1),C./cnorm));
    else
        C = max(variancelimit,C./cnorm);
    end

end

% get state 1 coefficients without intercept and return them in inverse
% filter form
A = [1 -B(1,(intercept+1):(p+intercept))];