%% Test Basic Functionality
rng(1994);
T = 100;
lambda = 5;
load('poissproc_hom_test.mat');
times = poissproc_hom(lambda, T);

assert(isequal(times, times_expected));