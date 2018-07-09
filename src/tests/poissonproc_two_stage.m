% ATTENTION: This script assumes that it is run from the src folder!
% Generate figures for a Poisson process with random starting time
%
% Copyright 2018 Lukas Koestler (TUM)

clear all;
close all;


%% plotting
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');

FontSize=12;
set(0,'defaultAxesFontSize',FontSize);
set(0,'defaultLegendFontSize',FontSize);
set(0,'defaultTextFontSize',FontSize);
set(0,'defaultLegendFontSize',FontSize);
set(0,'defaultLegendFontSize',FontSize);

%% parameters
num_samples = 1e6;
 
nu = 1000;
lambda= 1e-3;
sigma = 2;

lambda2 = lambda;
nu2 = nu;

% This change leads to E[N] = E[N2] (which is always true)
% but additionally Var(N) = Var(N2)
lambda2 = 4/3*lambda;
nu2 = 3/4*nu;
%% direct sampling
tic
N = poissproc_twostage_N(num_samples, nu, lambda, sigma);
t1 = toc;
fprintf('Time for direct sampling: %d sec per sample\n', t1/numel(N));


%% approximate sampling 
%N2 = zeros(1, num_samples);

tic
% for i = 1:numel(N2)
%     Nu = poissrnd(sigma*nu2);
%     N2(i) = numel(poissproc(@(t) Nu/sigma*lambda2*t, Nu*lambda2, sigma));
% end
N2 = poissrnd(sigma*nu2, 1, num_samples);
N2 = poissrnd(N2*lambda2*sigma/2);
t2 = toc;
fprintf('Time for approx sampling: %d sec per sample\n\n\n', t2/numel(N2));

%% Display
if lambda2 ~= lambda
    fprintf('--- Using changed nu, lambda for approx! ---\n');
end

fprintf('Mean (analytic, no approx): %d\n', lambda*nu*sigma^2 / 2);
fprintf('Mean (analytic, w/ approx): %d\n', lambda*nu*sigma^2 / 2);
fprintf('Mean (samples , no approx): %d\n', mean(N));
fprintf('Mean (samples , w/ approx): %d\n\n', mean(N2));

fprintf('Var  (analytic, no approx): %d\n', lambda*sigma^2*nu/4*(4/3*lambda*sigma + 2));
fprintf('Var  (analytic, w/ approx): %d\n', lambda2*sigma^2*nu2/4*(lambda2*sigma + 2));
fprintf('Var  (samples , no approx): %d\n', var(N));
fprintf('Var  (samples , w/ approx): %d\n\n\n', var(N2));

%% export table
if lambda2 ~= lambda
    fileID = fopen('../tex/paper/tables/poissproc02_changed.tex','w');
else
    fileID = fopen('../tex/paper/tables/poissproc02.tex','w');
end

fprintf(fileID,'Mean & No & Analytic & %d \\\\ \n', lambda*nu*sigma^2 / 2);
fprintf(fileID,'Mean & No & Samples & %d \\\\ \n', mean(N));
fprintf(fileID,'Mean & Yes & Analytic & %d \\\\ \n', lambda*nu*sigma^2 / 2);
fprintf(fileID,'Mean & Yes & Samples & %d \\\\ \n', mean(N2));

fprintf(fileID, '\\hline\\hline \n');
fprintf(fileID,'Variance & No & Analytic & %d \\\\ \n',...
    lambda*sigma^2*nu/4*(4/3*lambda*sigma + 2));
fprintf(fileID,'Variance & No & Samples & %d \\\\ \n', var(N));
fprintf(fileID,'Variance & Yes & Analytic & %d \\\\ \n',...
    lambda2*sigma^2*nu2/4*(lambda2*sigma + 2));
fprintf(fileID,'Variance & Yes & Samples & %d \\\\ \n', var(N2));

%% plot
tab1 = tabulate(N);
tab2 = tabulate(N2);

figure('Units', 'inches',...
    'position', [2,2,4,2.5],...
    'color', 'none',...
    'PaperPositionMode', 'auto');


bar(tab1(:,1), [tab1(:,3), tab2(:,3)]);
ylabel('Frequency in Percent');
xlabel('Number occurrences at $t=\sigma$'); 
legend({'Direct', 'Approximate'}, 'Color', 'none');
set(gca(), 'Color', 'none');

if lambda2 ~= lambda
    export_fig('-transparent', '../tex/paper/figures/poissproc02_changed_1.eps');
else
    export_fig('-transparent', '../tex/paper/figures/poissproc02_1.eps');
end

clf();
bar(tab1(:,1), [tab1(:,3), tab2(:,3)]);
ylabel('Frequency in Percent');
xlabel('Number occurrences at $t=\sigma$'); 
legend({'Direct', 'Approximate'}, 'Color', 'none');
set(gca(), 'YScale', 'log');
set(gca(), 'Color', 'none');
if lambda2 ~= lambda
    export_fig('-transparent', '../tex/paper/figures/poissproc02_changed_2.eps');
else
    export_fig('-transparent', '../tex/paper/figures/poissproc02_2.eps');
end
close all;
%% save workspace
if lambda2 ~= lambda
    save(sprintf('../data/poissonproc_two_stage_changed_1e%i.mat', round(log10(num_samples))));
else
    save(sprintf('../data/poissonproc_two_stage_1e%i.mat', round(log10(num_samples))));
end
%% Function Definitions

function N = poissproc_twostage_N(num_samples, nu, lambda, sigma)
% samples the number N for a two stage poisson process
% generates num_samples independent runs of the two-stage process

    N = zeros(1, num_samples);
    for i = 1:num_samples
        us = poissproc_hom(nu, sigma);
        N(i) = sum(poissrnd((sigma-us)*lambda));
    end
end

