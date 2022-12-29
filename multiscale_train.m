function B = multiscale_train(X,lags)
% B = multiscale_train(X,lags)
% Learn a multi-scale autoregressive (AR) filter to represent the time evolution
% of the columns of matrix X.
%
% Inputs:
% X - A vector time series, with observation vectors as columns.
% lags - A cell array whose each element is a vector containing the lag
% structure of one AR filter, the length of the cell array determines the
% number of filters.
%
% Outputs:
% B - Filter coefficients estimated using least squares. Filter
% coefficients are in the form of sparse vector autoregression matrices.
%
% Example:
% X = randn(39,10000);
% B = multiscale_train(X,{1:10,2:2:20,3:3:30,4:4:40,5:5:50});
% Y = randn(39,1000);
% [Yh,e] = multiscale_filter(Y,B,{1:10,2:2:20,3:3:30,4:4:40,5:5:50});
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

persistent buffer_function_exists;

% check if the "buffer" function is implemented on this platform
if isempty(buffer_function_exists)
    buffer_function_exists = exist('buffer','file');
end

for i1=1:length(lags) % for each filter
    p = length(lags{i1}); % order of filter
    d = lags{i1}; % lag structure of filter
    if isscalar(d)
        d = (1:p)*d;
    end

    [K T] = size(X); % K=number of variables, T=number of frames

    % an order p vector autoregressive filter for K-dimensional vector sequence 
    % with intercept vector as the first column
    B{i1} = zeros(K,(K*p+1));

    A = zeros(K,p+1);
    for i2=1:K % estimate autoregressive model for each variable
        if buffer_function_exists
            D = flipud(buffer([0 X(i2,1:end-1)],d(p),d(p)-1)); % d(p) is maximum lag
        else
            D = zeros(d(p),T);
            for i3=1:d(p)
                D(i3,(i3+1):T) = X(i2,1:(T-i3));
            end
        end
        Z = [ones(1,T);D(d,:)]; % concatenated with ones because of the intercept term
        A(i2,:) = (Z*Z')\Z*X(i2,:)'; % least squares
    end
    % convert the estimated AR models to vector autoregressive (VAR) form
    % (for simple and efficient computation in Matlab)
    % the VAR matrices will be sparse because all the prediction
    % coefficients related to other variables than the predicted one will
    % be zero
    B{i1}(:,1) = A(:,1); % first column: intercept terms
    for i2=1:p
        % the columns 2...K*p+1 contain lagged AR coefficients
        B{i1}(:,((i2-1)*K+2):(i2*K+1)) = diag(A(:,i2+1));
    end
end