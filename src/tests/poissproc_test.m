%% Test Basic Functionality
rng(1994);
T = 100;
lambda = @(t) t;
lambda_max = T;

load('poissproc_test.mat');
times = poissproc(lambda, lambda_max, T);

assert(isequal(times, times_expected));


%% Test fallback to poissproc_hom
rng(1994);
T = 100;
lambda = 5;
lambda_max = 5;

times1 = poissproc(lambda, lambda_max, T);
rng(1994);
times2 = poissproc_hom(lambda, T);

assert(isequal(times1, times2));