function [N, M, ind, times] = poissproc2_slow(eta, lambda, sigma, num_samples)
%POISSPROC2_SLOW Two stage Poisson process slow simulation
%
% Two stage Poisson process approximate sampling
%
%   [N, times, ind] = POISSPROC2_SLOW(eta, lambda, sigma, num_samples)
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
%       ind: starting index of output times. times(ind(i):ind(i+1)-1)
%           belongs to sample i. This is unhandy but good for vectorization
%       times: arrival times of the second stage process for all samples,
%           not sorted, use ind to sort
%
%   DESCRIPTION:
%       'Direct' Computation. Much slower than poissproc2 but 'without'
%       approximations
%       The output format is very unintuitive but leads to speed up due to
%       vectorization
%
%
% Copyright 2018 Lukas Koestler (TUM)

%% preliminary work
etainv = inv(eta/eta(sigma));
m = cumsum(lambda);
minv = inv(m);

% expected number of occurences for the stage-one process
EN = eta(sigma);

%% First stage
N = poissrnd(EN, 1, num_samples);
indTs = repelem(1:num_samples, N);
Ts = etainv(rand(1, sum(N)));

Ms = poissrnd(m(sigma-Ts));
ind_M = find(Ms > 0);
us = m(sigma - Ts(ind_M)) .* rand(1, numel(ind_M));

times = minv(us)+Ts(ind_M);
indTs2 = indTs(ind_M);

tab = tabulate(indTs2);
M = [tab(:,2).', zeros(1, num_samples-tab(end,1))];
ind = zeros(1, num_samples+1);
ind(2:end) = cumsum(M)+1;
ind(1)=1;

clearvars indTs Ts Ms ind_M indTs2 tab