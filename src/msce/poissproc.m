 function times = poissproc(lambda, lambda_max, T)
%POISSPROC Poisson process arrival times
%
% Generate the arrival times of a general poisson process
%
%   times = POISSPROC(lambda, lambda_max, T)
%
%   INPUT:
%       lambda: rate of process > 0, function or scalar
%               if lambda is a function it must be vectorized
%       lambda_max: maximum of the rate, scalar
%       T: end time of process > 0, scalar
%   OUTPUT:
%       times: arrival times. >= 0, < T, vector
%
%   DESCRIPTION:
%       Uses thinning. If lambda is a scalar, then a faster algorithm
%       (poissproc_hom) is used.
%
% See also poissproc_hom
%
% Copyright 2018 Lukas Koestler (TUM)

if isnumeric(lambda) && isscalar(lambda)
    times = poissproc_hom(lambda, T);
else
    times = poissproc_hom(lambda_max, T);
    
    % filter
    u = rand(size(times));
    lambda_at_times = lambda(times);
    if max(lambda_at_times) > lambda_max
        error('lambda(t) <= lambda_max does not hold!');
    end
    ind_keep = (u*lambda_max) <= lambda_at_times;
    times = times(ind_keep);
end
end


