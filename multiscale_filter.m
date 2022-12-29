function [Xh,res] = multiscale_filter(X,B,lags)
% [Xh,res] = multiscale_filter(X,B,lags)
%
% Inputs:
% X - A sequence of column vectors.
% B - A cell array of vector autoregressive parameter matrices (one for each filter).
% lags - A cell array whose each element is a vector containing the lag
% structure of one AR filter. The length of this cell array should be the same as (or not less
% than) cell array B.
%
% Outputs:
% Xh - Multi-scale autoregressive filtered version of X.
% res - Prediction residual matrix, same size as X.
%
% The following paper should be cited in case of this method being
% used in a publication:
% J. Pohjalainen and P. Alku, "Multi-scale modulation filtering in
% automatic detection of emotions in telephone speech", in Proc. 39th
% International Conference on Acoustics, Speech, and Signal Processing
% (ICASSP), Florence, Italy, May 4-9 2014.
% Please also see that paper for a detailed description of the method and
% its possible usage. The paper is available for download at
% http://users.tkk.fi/~jpohjala/.
%
% Author: Jouni Pohjalainen, Department of Signal Processing and Acoustics,
% Aalto University, 2014.
% Mail: jpohjala@acoustics.hut.fi
% The algorithm can be freely used for research purposes. Please send your
% bug reports or other comments to the author.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[K,T] = size(X);

for i1=1:length(B)
    % order p is determined based on the dimensions of the vector autoregressive matrix
    p = (size(B{i1},2)-1)/K;
    d = lags{i1};

    % matrix Z contains all the variables so that we can
    % directly multiply it by the vector autoregression matrix
    Z = zeros(K*p+1,T);
    Z(1,:) = ones(1,T); % the first row is ones because of the intercept terms
    for i2=1:p
        Z(((i2-1)*K+2):(1+i2*K),:) = [zeros(K,min([d(i2) T])) X(:,1:end-d(i2))];
    end

    res(:,:,i1) = X-B{i1}*Z; % residual matrix for this filter
    if i1==length(B) % if all filters have been used, choose the best one for each time instant and variable
        [cols,rows] = meshgrid(1:size(res,2),1:size(res,1));
        [tmp,idx] = min(res.^2,[],3); % for each prediction, choose residual whose squared error is minimum
        res = res(sub2ind(size(res),rows,cols,idx));
    end
end

Xh = X - res; % get final prediction by subtracting the minimum residual from original data