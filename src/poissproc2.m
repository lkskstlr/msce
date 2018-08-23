 function [N, M, ind, times] = poissproc2(eta, lambda, sigma, num_samples)
%POISSPROC2 Two stage Poisson process
%
% Two stage Poisson process approximate sampling
%
%   [N, times, ind] = poissproc2(eta, lambda, sigma, num_samples)
%
%   INPUT:
%       eta: mean value function of the first process, chebfun,
%       monotonically increasing, the domain MUST be [0, sigma]
%       lambda: rate of the second process > 0, chebfun,
%               integral(lambda, 0, sigma) << 1 or the approximation is
%               invalid, the domain MUST be [0, sigma]
%       sigma: 0 <= t <= sigma, interval length
%       num_samples: number of samples of the two-stage poisson process
%           to simulate. Because of MATLABs vectorized nature it is
%           advisable to simulate all runs at once
%   OUTPUT:
%       N: number of occurences in the first stage for all samples,
%           vector with length num_samples
%       M: number of occurences in the second stage for all samples,
%           vector with length num_samples
%       ind: starting index of output times. times(ind(i):ind(i+1)-1)
%           belongs to sample i. This is unhandy but good for vectorization
%       times: arrival times of the second stage process for all samples,
%           not sorted, use ind to sort
%
%   DESCRIPTION:
%       See Theorem 1.3
%       The output format is very unintuitive but leads to speed up due to
%       vectorization
%
%
% Copyright 2018 Lukas Koestler (TUM)


%% preliminary work

% expected number of occurences for the stage-one process
EN = eta(sigma);

% mean value function of mu
h = conv(eta, lambda, 'same');

% normalize
hmax = h(sigma);
h = h/hmax;

%% sample number of occurences
% This can be heavily improved by
% N = EN + round(sqrt(EN)*randn(1, num_samples));
N = poissrnd(EN, 1, num_samples);
M = poissrnd(N*hmax/EN);

ind = zeros(1, num_samples+1);
ind(2:end) = cumsum(M)+1;
ind(1)=1;

%% sample arrival times
times = helper_sample_h(h, sigma, 0.001*sigma, sum(M));
end


