% ATTENTION: This script assumes that it is run from the src folder!
% Generate figures for a Poisson process with random starting time
%
% Copyright 2018 Lukas Koestler (TUM)

clear all;
close all;


%% define rates
nu = chebfun('1000*(1+0.5*sin(4*pi*x))', [0, 1]); %first rate
eta = cumsum(nu); %firts mean value function
lambda = chebfun('0.001*(1-0.02*x^2)', [0, 1]); %second rate
sigma = 1; %interval length


%% simulation without approximation
num_samples = 1e4;
tic;
[N, M, ind, times] = poissproc2_slow(eta, lambda, sigma, num_samples);
fprintf('No   Approximation: %d sec/sample\n',toc/num_samples);

%% simulation with approximation
num_samples = 1e6;
tic;
[N2, M2, ind2, times2] = poissproc2(eta, lambda, sigma, num_samples);
fprintf('With Approximation: %d sec/sample\n',toc/num_samples);

%% plot
figure();
title('Arrival times after two stage process');
histogram(times2, 'Normalization', 'pdf');
hold on
histogram(times, 'Normalization', 'pdf');
legend('With Approximation', 'No Approximation');