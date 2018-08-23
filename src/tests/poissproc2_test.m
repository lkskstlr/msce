%% Test Basic Functionality
sigma = 50;
eta = cumsum(chebfun(@(t) 1.0 + 0.2*sin(t), [0, sigma]));
lambda = chebfun(@(t) (1 + 0.2*cos(t))/(sigma*15), [0, sigma]);
num_samples = 1e3;

[N, M, ind, times] =...
    poissproc2(eta, lambda, sigma, num_samples);
[N_slow, M_slow, ind_slow, times_slow] =...
    poissproc2_slow(eta, lambda, sigma, num_samples);
