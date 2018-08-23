clear all;
close all;

%% Parameters
N = 1e5;
sigma = 50;
model = crc_mkmodel();


%% simulate to pre-screen
tic;
[N2, N3, N4, N_polyp] = crc_simulation(model, sigma, N);
time_sim = toc;
fprintf('Simulation: %d sec/sample, %d samples/sec \n',...
    time_sim/N,...
    N/time_sim);

%% global plot settings
cols_ = cbrewer('qual', 'Set1', 5);
colb = cols_(3,:);
col1 = cols_(4,:);
col2 = cols_(2,:);
col3 = cols_(1,:);

col5 = cols_(5,:);

col_back = [35, 55, 59]/255;
col_w = [250, 250, 250]/255;

set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');

FontSize=12;
set(0,'defaultAxesFontSize',FontSize);
set(0,'defaultLegendFontSize',FontSize);
set(0,'defaultTextFontSize',FontSize);
set(0,'defaultLegendFontSize',FontSize);
set(0,'defaultLegendFontSize',FontSize);

%% Screening
thres1 = 1e3;
thres2 = 1e4;
thres3 = 1e5;
screening1 = zeros(1, N);
screening2 = zeros(1, N);
screening3 = zeros(1, N);

for i = 1:N
    screening1(i) = any(N_polyp{i} > thres1);
    screening2(i) = any(N_polyp{i} > thres2);
    screening3(i) = any(N_polyp{i} > thres3);
end

%% hazard functions
xx = linspace(50, 100, 1e3+1);
i80 = find(xx == 80);

LogS = zeros(N, length(xx));
H = zeros(N, length(xx));

logS1 = crc_mklogsurvival1(model, 100-sigma);
logS2 = crc_mklogsurvival2(model, 100-sigma);
logS3 = crc_mklogsurvival3(model, 100-sigma);
logS4 = crc_mklogsurvival4(model, 100-sigma);

H1 = msce_hazard1(model, xx-sigma);
H2 = msce_hazard2(model, xx-sigma);
H3 = msce_hazard3(model, xx-sigma);
H4 = msce_hazard4(model, xx-sigma);

LogS1 = logS1(xx-sigma);
LogS2 = logS2(xx-sigma);
LogS3 = logS3(xx-sigma);
LogS4 = logS4(xx-sigma);

for i = 1:N
    LogS(i,:) = model.X*LogS4 + N2(i)*LogS3 + N3(i)*LogS2 + N4(i)*LogS1;
    H(i,:) = model.X*H4 + N2(i)*H3 + N3(i)*H2 + N4(i)*H1;
end

Sneg = mean(exp(LogS(screening1==0, :)), 1);
Sall = mean(exp(LogS), 1);

Hall = sum( exp(LogS) .* H, 1) ./ sum(exp(LogS));
Hneg1 = sum( exp(LogS(screening1==0,:)) .* H(screening1==0, :), 1) ./ sum( exp(LogS(screening1==0,:)), 1);
Hneg2 = sum( exp(LogS(screening2==0,:)) .* H(screening2==0, :), 1) ./ sum( exp(LogS(screening2==0,:)), 1);
Hneg3 = sum( exp(LogS(screening3==0,:)) .* H(screening3==0, :), 1) ./ sum( exp(LogS(screening3==0,:)), 1);

Hpos2 = sum( exp(LogS(screening2==1,:)) .* H(screening2==1, :), 1) ./ sum( exp(LogS(screening2==1,:)), 1);

%% Intervention complete
%thres 1
LogSinter = zeros(N, length(xx));
Hinter = zeros(N, length(xx));
thres = thres1;

j = 1;
for i = 1:N
    ind = N_polyp{i} <= thres;
    if any(ind == 0)
        N3_loc = sum(ind);
        N4_loc = sum(N_polyp{i}(ind));
        LogSinter(j,:) = model.X*LogS4 + N2(i)*LogS3 + N3_loc*LogS2 + N4_loc*LogS1;
        Hinter(j,:) = model.X*H4 + N2(i)*H3 + N3_loc*H2 + N4_loc*H1;
        j=j+1;
    end
end

LogSinter1 = LogSinter(1:j-1, :);
Hinter1 = Hinter(1:j-1, :);

% thres 2
LogSinter = zeros(N, length(xx));
Hinter = zeros(N, length(xx));
thres = thres2;

j = 1;
for i = 1:N
    ind = N_polyp{i} <= thres;
    if any(ind == 0)
        N3_loc = sum(ind);
        N4_loc = sum(N_polyp{i}(ind));
        LogSinter(j,:) = model.X*LogS4 + N2(i)*LogS3 + N3_loc*LogS2 + N4_loc*LogS1;
        Hinter(j,:) = model.X*H4 + N2(i)*H3 + N3_loc*H2 + N4_loc*H1;
        j=j+1;
    end
end

LogSinter2 = LogSinter(1:j-1, :);
Sinter2 = mean(exp(LogSinter2), 1);
Hinter2 = Hinter(1:j-1, :);

% thres 3
LogSinter = zeros(N, length(xx));
Hinter = zeros(N, length(xx));
thres = thres3;

j = 1;
for i = 1:N
    ind = N_polyp{i} <= thres;
    if any(ind == 0)
        N3_loc = sum(ind);
        N4_loc = sum(N_polyp{i}(ind));
        LogSinter(j,:) = model.X*LogS4 + N2(i)*LogS3 + N3_loc*LogS2 + N4_loc*LogS1;
        Hinter(j,:) = model.X*H4 + N2(i)*H3 + N3_loc*H2 + N4_loc*H1;
        j=j+1;
    end
end

LogSinter3 = LogSinter(1:j-1, :);
Hinter3 = Hinter(1:j-1, :);

HHinter1 = sum( exp(LogSinter1) .* Hinter1, 1) ./ sum(exp(LogSinter1), 1);
HHinter2 = sum( exp(LogSinter2) .* Hinter2, 1) ./ sum(exp(LogSinter2), 1);
HHinter3 = sum( exp(LogSinter3) .* Hinter3, 1) ./ sum(exp(LogSinter3), 1);


%% Intervention without progenitor
%thres 1
LogSinterin = zeros(N, length(xx));
Hinterin = zeros(N, length(xx));
thres = thres1;

j = 1;
for i = 1:N
    ind = N_polyp{i} <= thres;
    if any(ind == 0)
        N4_loc = sum(N_polyp{i}(ind));
        LogSinterin(j,:) = model.X*LogS4 + N2(i)*LogS3 + N3(i)*LogS2 + N4_loc*LogS1;
        Hinterin(j,:) = model.X*H4 + N2(i)*H3 + N3(i)*H2 + N4_loc*H1;
        j=j+1;
    end
end

LogSinterIn1 = LogSinterin(1:j-1, :);
HinterIn1 = Hinterin(1:j-1, :);

% thres 2
LogSinterin = zeros(N, length(xx));
Hinterin = zeros(N, length(xx));
thres = thres2;

j = 1;
for i = 1:N
    ind = N_polyp{i} <= thres;
    if any(ind == 0)
        N4_loc = sum(N_polyp{i}(ind));
        LogSinterin(j,:) = model.X*LogS4 + N2(i)*LogS3 + N3(i)*LogS2 + N4_loc*LogS1;
        Hinterin(j,:) = model.X*H4 + N2(i)*H3 + N3(i)*H2 + N4_loc*H1;
        j=j+1;
    end
end

LogSinterIn2 = LogSinterin(1:j-1, :);
HinterIn2 = Hinterin(1:j-1, :);

% thres 3
LogSinterin = zeros(N, length(xx));
Hinterin = zeros(N, length(xx));
thres = thres3;

j = 1;
for i = 1:N
    ind = N_polyp{i} <= thres;
    if any(ind == 0)
        N4_loc = sum(N_polyp{i}(ind));
        LogSinterin(j,:) = model.X*LogS4 + N2(i)*LogS3 + N3(i)*LogS2 + N4_loc*LogS1;
        Hinterin(j,:) = model.X*H4 + N2(i)*H3 + N3(i)*H2 + N4_loc*H1;
        j=j+1;
    end
end

LogSinterIn3 = LogSinterin(1:j-1, :);
HinterIn3 = Hinterin(1:j-1, :);

HHinterIn1 = sum( exp(LogSinterIn1) .* HinterIn1, 1) ./ sum(exp(LogSinterIn1), 1);
HHinterIn2 = sum( exp(LogSinterIn2) .* HinterIn2, 1) ./ sum(exp(LogSinterIn2), 1);
HHinterIn3 = sum( exp(LogSinterIn3) .* HinterIn3, 1) ./ sum(exp(LogSinterIn3), 1);

%% Lifetime risk with intervention
LogSinterLife = zeros(N, length(xx));
HinterLife = zeros(N, length(xx));
thres = thres2;

j = 1;
for i = 1:N
    ind = N_polyp{i} <= thres;
    if any(ind == 0)
        N4_loc = sum(N_polyp{i}(ind));
        N4_loc = N4_loc + 0.1*thres*sum(ind==0);
        LogSinterLife(j,:) = model.X*LogS4 + N2(i)*LogS3 + N3(i)*LogS2 + N4_loc*LogS1;
        HinterLife(j,:) = model.X*H4 + N2(i)*H3 + N3(i)*H2 + N4_loc*H1;
        j=j+1;
    end
end

LogSinterLife = LogSinterLife(1:j-1, :);
HinterLife = HinterLife(1:j-1, :);

SinterLife = mean(exp(LogSinterLife), 1);
%% Table progenitor cells
tab = tabulate(N3);
fileID = fopen('tab_apc_mm.tex','w');
for i = 1:size(tab, 1)
    fprintf(fileID, '%i & %i & %.2f \\%% \\\\ \n', tab(i,1), tab(i,2), tab(i,3));
end

%% size distribution polyps
fig = figure('Color', 'none',...
    'Units', 'inches', 'Position', [5, 5, 4, 2.5],...
    'PaperPositionMode', 'auto');
edges = linspace(0,7,15);
h_hist = histogram(log10(cat(2,N_polyp{:})), edges,...
    'Normalization', 'probability');
h_hist.FaceAlpha = 1;
h_hist.LineWidth = 1.5;
h_hist.EdgeColor = col_w;%'none';
h_hist.FaceColor = col_back;
hold on
line(log10([1e5 1e5]), ylim(), 'Color', col3, 'LineWidth', 2);
line(log10([1e4 1e4]), ylim(), 'Color', col2, 'LineWidth', 2);
line(log10([1e3 1e3]), ylim(), 'Color', col1, 'LineWidth', 2);

set(gca(), 'Color', 'none');

legend({'Size', '$10^5$', '$10^4$', '$10^3$'}, ...
     'Location', 'NorthWest');

xticks((0:6));
xticklabels_ = cell(1,7);
for i = 0:6
    xticklabels_{i+1} = sprintf('10^%i', i);
end
xticklabels(xticklabels_);
xlabel('Number of cells in polyp');
ylabel('Fraction');

% print -depsc2 test.eps
export_fig -transparent polyp_size_dist.eps

close(fig)

%% Mortality
fileID = fopen('tab_lifetime_risk.tex','w');
fprintf(fileID,'Background & & %.2f \\%% \\pause \\\\ \n', 100*(1 - Sall(i80)));

fprintf(fileID, '\\hline \n');
fprintf(fileID,'Neg. Screen & $10^5$ & %.2f \\%% \\\\ \n',...
    100*(1 - mean(exp(LogS(screening3==0,i80)))));
fprintf(fileID,'Neg. Screen & $10^4$ & %.2f \\%% \\\\ \n',...
    100*(1 - mean(exp(LogS(screening2==0,i80)))));
fprintf(fileID,'Neg. Screen & $10^3$ & %.2f \\%% \\pause \\\\ \n',...
    100*(1 - mean(exp(LogS(screening1==0,i80)))));

fprintf(fileID, '\\hline \n');
fprintf(fileID,'Pos. Screen & $10^5$ & %.2f \\%% \\\\ \n',...
    100*(1 - mean(exp(LogS(screening3==1,i80)))));
fprintf(fileID,'Pos. Screen & $10^4$ & %.2f \\%% \\\\ \n',...
    100*(1 - mean(exp(LogS(screening2==1,i80)))));
fprintf(fileID,'Pos. Screen & $10^3$ & %.2f \\%% \\pause \\\\ \n',...
    100*(1 - mean(exp(LogS(screening1==1,i80)))));

fprintf(fileID, '\\hline \n');
fprintf(fileID,'Realistic Intervention & $10^4 \\rightarrow 10^3$ & %.2f \\%% \\\\ \n',...
    100*(1 - SinterLife(i80)));
fprintf(fileID,'Complete Intervention & $10^4 \\rightarrow \\phantom{1}0$ & %.2f \\%% \\\\ \n',...
    100*(1 - Sinter2(i80)));
%% Figure 3 Paper
fig = figure('Color', 'none',...
    'Units', 'inches', 'Position', [5, 5, 4, 2.5],...
    'PaperPositionMode', 'auto');
semilogy(xx, Hall, 'Color', colb)
hold on
semilogy(xx, Hneg3, 'Color', col3)
semilogy(xx, Hneg2, '--', 'Color', col2)
semilogy(xx, Hneg1, 'Color', col1)
ylim([1e-6, 1e-1]);
yticks(10.^(-7:-1));

set(gca(), 'Color', 'none');

legend({'background', 'threshold $10^5$', 'threshold $10^4$', 'threshold $10^3$'}, ...
     'Location', 'SouthEast');


xlabel('Age');
ylabel('Hazard');

% print -depsc2 test.eps
export_fig -transparent paper_figure3.eps

close(fig)


%% Figure 4 Paper
fig = figure('Color', 'none',...
    'Units', 'inches', 'Position', [5, 5, 4, 2.5],...
    'PaperPositionMode', 'auto');
semilogy(xx, Hall, 'Color', colb)
hold on
semilogy(xx, Hpos2, 'Color', col5);
semilogy(xx, Hneg2, 'Color', col2);
ylim([1e-5, 1e+1]);
yticks(10.^(-6:+2));
set(gca(), 'Color', 'none');

legend({'background', 'Pos. Screen ($10^4$)', 'Neg. Screen ($10^4$)'}, ...
     'Location', 'SouthEast');


xlabel('Age');
ylabel('Hazard');

% print -depsc2 test.eps
export_fig -transparent paper_figure4.eps

close(fig)

%% Figure 5 Paper
fig = figure('Color', 'none',...
    'Units', 'inches', 'Position', [5, 5, 4, 2.5],...
    'PaperPositionMode', 'auto');
semilogy(xx, Hall, 'Color', colb)
hold on
semilogy(xx, HHinter3, 'Color', col3);
semilogy(xx, HHinter2, '--', 'Color', col2);
semilogy(xx, HHinter1, 'Color', col1);
ylim([1e-6, 1e-1]);
yticks(10.^(-7:-1));

set(gca(), 'Color', 'none');

legend({'background', 'threshold $10^5$', 'threshold $10^4$', 'threshold $10^3$'}, ...
     'Location', 'SouthEast');


xlabel('Age');
ylabel('Hazard');

% print -depsc2 test.eps
export_fig -transparent paper_figure5.eps

close(fig)


%% Figure 6 Paper
fig = figure('Color', 'none',...
    'Units', 'inches', 'Position', [5, 5, 4, 2.5],...
    'PaperPositionMode', 'auto');
semilogy(xx, Hall, 'Color', colb)
hold on
semilogy(xx, HHinterIn3, 'Color', col3);
semilogy(xx, HHinterIn2, '--', 'Color', col2);
semilogy(xx, HHinterIn1, 'Color', col1);
ylim([1e-6, 1e-1]);
yticks(10.^(-7:-1));

set(gca(), 'Color', 'none');

legend({'background', 'threshold $10^5$', 'threshold $10^4$', 'threshold $10^3$'}, ...
     'Location', 'SouthEast');


xlabel('Age');
ylabel('Hazard');

% print -depsc2 test.eps
export_fig -transparent paper_figure6.eps

close(fig)








