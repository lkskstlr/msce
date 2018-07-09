% ATTENTION: This script assumes that it is run from the src folder!
% Generate figures for a Poisson process with random starting time
%
% Copyright 2018 Lukas Koestler (TUM) 


clear variables;
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
num_samples = 1e7;

sigma = 1;
p_rnd = @(n) rand(1,n);


figure('Units', 'inches',...
    'position', [2,2,4,2.5],...
    'color', 'none',...
    'PaperPositionMode', 'auto');
for j = 0:2
    lambda = 10^(-j);
    N = poissproc_randtau_N(num_samples, lambda, sigma, p_rnd);
    tab = tabulate(N);

    clf();
    bar(tab(:,1), tab(:,3)/100);
    hold on;
    plot(tab(:,1), poisspdf(tab(:,1), lambda*sigma^2/2), '+r-');
    set(gca(), 'Color', 'none');
    export_fig('-transparent', sprintf('../tex/paper/figures/poissproc01_%i_%i.eps', j, 1));

    clf();
    bar(tab(:,1), tab(:,3)/100);
    hold on;
    plot(tab(:,1), poisspdf(tab(:,1), lambda*sigma^2/2), '+r-');
    set(gca(), 'YScale', 'log');
    set(gca(), 'Color', 'none');
    export_fig('-transparent', sprintf('../tex/paper/figures/poissproc01_%i_%i.eps', j, 2));

end

close all;

%% Function Definitions

function N = poissproc_randtau_N(num_samples, lambda, sigma, p_rnd)
% samples the number N for a poisson process with random starting time
    tau = p_rnd(num_samples);
    N = poissrnd(lambda*(sigma-tau));
end

