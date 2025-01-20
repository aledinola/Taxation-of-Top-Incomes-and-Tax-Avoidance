%% Show comparative statics results
%% #VC# V12
% Last updated by Alessandro on 2021-04-04

clear; clc; close all
format long g

% - Directory where figures are saved
FigDir = fullfile('results','figures'); 

% - Flag 0/1 if want to save plots as png
do_save = 1;

% - Attach suffix to each figure for experiment, or leave it as ''
exp1 = '_small_grid';
exp2 = '_large_grid';

% Create folder for figures if not existent
if ~exist(FigDir, 'dir')
    mkdir(FigDir)
end

%% Load data from txt/excel/manual copying

% Experiment 1: nl=20
% Experiment 2: nl=100

% Add 1.41 for plots

results(1).lambda = [1.41
1.49
1.5
1.500001
1.505
1.51
1.52
1.55
1.6
1.7
1.8
1.9
2];


results(1).share_ec = [0
0
0
0.053895
0.114472
0.11429
0.114093
0.115123
0.114008
0.121216
0.162515
0.405477
0.626984];

results(1).share_entre = [0.127774
0.127774
0.127774
0.123891
0.124857
0.124987
0.125075
0.125096
0.125526
0.125306
0.125925
0.125536
0.125886];

% Now nl is 100
results(2).lambda = [1.41
1.49
1.5
1.500001
1.51
1.52
1.55
1.8
1.9
1.925
1.95
2
2.05];



results(2).share_ec = [0
0
0
0
0
0
0
0
0.028875
0.583205
0.584418
0.608369
0.61554];



results(2).share_entre = [0.12693
0.126941
0.126941
0.126941
0.12691
0.126833
0.126852
0.128613
0.129013
0.125884
0.125701
0.125763
0.126265];


%% Plots for nl = 20
ind = 0;


ind = ind+1;
figure(ind)
plot(results(1).lambda,results(1).share_ec,'-o','linewidth',2)
hold on
plot(results(1).lambda,results(1).share_entre,'-o','linewidth',2)
legend('Share of C | Entre','Total share of Entre','Location','best','fontsize',14)
xlabel('\lambda_{C}','fontsize',14)
ylabel('Share','fontsize',14)
title('Share of C-corps out of Entrepr., NL = 20')
hold off
if do_save==1; print([fullfile(FigDir,'lambda_share_ec'),exp1],'-dpng'); end

ind = ind+1;
figure(ind)
subplot(1,2,1)
plot(results(1).lambda,results(1).share_ec,'-o','linewidth',2)
legend('Share of C | Entre','Location','best','fontsize',14)
xlabel('\lambda_{C}','fontsize',14)
ylabel('Share','fontsize',14)
title('Share of C-corps out of Entrepr., NL=100')

subplot(1,2,2)
plot(results(1).lambda,results(1).share_entre,'-o','linewidth',2)
ylim([0, 0.7])
legend('Total share of Entre','Location','best','fontsize',14)
xlabel('\lambda_{C}','fontsize',14)
ylabel('Share','fontsize',14)
title('Share of Entrepr., NL = 20')
if do_save==1; print([fullfile(FigDir,'sub_lambda_share_ec'),exp1],'-dpng'); end

%% Plots for NL = 100
ind = ind+1;
figure(ind)
plot(results(2).lambda,results(2).share_ec,'-o','linewidth',2)
hold on
plot(results(2).lambda,results(2).share_entre,'-o','linewidth',2)
legend('Share of C | Entre','Total share of Entre','Location','best','fontsize',14)
xlabel('\lambda_{C}','fontsize',14)
ylabel('Share','fontsize',14)
title('Share of C-corps out of Entrepr., NL=100')
hold off
if do_save==1; print([fullfile(FigDir,'lambda_share_ec'),exp2],'-dpng'); end

ind = ind+1;
figure(ind)
subplot(1,2,1)
plot(results(2).lambda,results(2).share_ec,'-o','linewidth',2)
legend('Share of C | Entre','Location','best','fontsize',14)
xlabel('\lambda_{C}','fontsize',14)
ylabel('Share','fontsize',14)
title('Share of C-corps out of Entrepr., NL=100')

subplot(1,2,2)
plot(results(2).lambda,results(2).share_entre,'-o','linewidth',2)
ylim([0, 0.7])
legend('Total share of Entre','Location','best','fontsize',14)
xlabel('\lambda_{C}','fontsize',14)
ylabel('Share','fontsize',14)
title('Share of Entrepr., NL=100')
if do_save==1; print([fullfile(FigDir,'sub_lambda_share_ec'),exp2],'-dpng'); end



