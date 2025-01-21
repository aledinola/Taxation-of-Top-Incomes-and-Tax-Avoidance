function fun_plot_compstat(base, CF2, file_root,do_save)
% PURPOSE
% fun_plot_compstat makes plots of results A in files with stem filestem
% INPUTS
% A         : struct
% filestem  : character with root of the filename
% do_save   : flag 0/1

if ~isstruct(base)
    error('Input argument base must be a structure')
end
if ~isstruct(CF2)
    error('Input argument CF2 must be a structure')
end
if ~ischar(file_root)
    error('Input argument file_root must be a character variable - string')
end

% Plotting options
mylw     = 4;  % line width
myfs   = 18; % font size
mybw = 1;
mybaredge = 3;
myfc1 = [0.03 0.05 0.29]; % facecolor for bar plots
myfc2 = [ 0.5    0.75    0.98];
myfc3 = [0.89 0.92 0.97];
mylc1 = [0 0 0];
mylc2 = [0, 0.4470, 0.7410];
mylc3 = [0.8500, 0.3250, 0.0980];
mylp1 = '-';
mylp2 = '--';
mylp3 = '-.';


% tau_h grid
tau_h_grid = [1:38];


% Baseline tax rate (0.396)
top_tax = 0.396;
base_toptax_ind = find(base.tau_h>=top_tax,1);
CF2_toptax_ind = find(CF2.tau_h>=top_tax,1);

% Total taxes paid by young pop
base.taxes_young = sum(base.taxes_occ,2);
CF2.taxes_young = sum(CF2.taxes_occ,2);

%% Calculate elasticity of occupational choice wrt tau_h
% percentage change in tau_h
pct_tau_h_base = 100*(base.tau_h(2)-base.tau_h(1))/base.tau_h(1);
pct_tau_h_CF2 = 100*(CF2.tau_h(2)-CF2.tau_h(1))/CF2.tau_h(1);

% percentage change in share of entre.
pct_share_entre_base = 100*(base.share_entre_act(2)-base.share_entre_act(1))/base.share_entre_act(1);
pct_share_entre_CF2 = 100*(CF2.share_entre_act(2)-CF2.share_entre_act(1))/CF2.share_entre_act(1);

% elasticity
elas_base = pct_share_entre_base/pct_tau_h_base;
elas_CF2 = pct_share_entre_CF2/pct_tau_h_CF2;

disp([elas_base elas_CF2])

% elas_base = 0.0146092996030178
% elas_CF2 = 0.116985270403935


%% Make plot for SUERF figure 3
myfs = myfs - 4;
figure('Position', [1, 1, 750, 500])
t=tiledlayout(2,3);
t.TileSpacing = 'compact';
nexttile
%share of C-corporations
plot(base.tau_h(tau_h_grid),100*base.share_EC_entre(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.share_EC_entre(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[-1,30];
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
ylabel('Share of entrepreneurs (%)','fontsize',myfs)
title(["(a) Share of C-corp.",newline],'FontSize',myfs)

% aggregate capital
nexttile
plot(base.tau_h(tau_h_grid),(base.K(tau_h_grid)/base.K(base_toptax_ind)),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),(CF2.K(tau_h_grid)/CF2.K(CF2_toptax_ind)),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim = [0.85,1.2];
% xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
ylabel('1 = Baseline \tau_h')
title(["(b) Agg. capital",newline],'FontSize',myfs)

% aggregate output
nexttile
plot(base.tau_h(tau_h_grid),(base.Y(tau_h_grid)/base.Y(base_toptax_ind)),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),(CF2.Y(tau_h_grid)/CF2.Y(CF2_toptax_ind)),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim = [0.85,1.2];
xlabel('\tau_h','fontsize',myfs)
ylabel('1 = Baseline \tau_h')
title(["(c) Agg. output",newline],'FontSize',myfs)

% top 1 % income share
nexttile
name = 'inc_share_top1';
plot(base.tau_h(tau_h_grid),100*base.(name)(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.(name)(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[12,18];
ylabel('Share of total income (%)','fontsize',myfs)
xlabel('\tau_h','fontsize',myfs)
title(["(d) Top 1% income share",newline],'FontSize',myfs)

% top 1% wealth share
nexttile
name = 'wealth_share_top1';
plot(base.tau_h(tau_h_grid),100*base.(name)(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.(name)(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[17,23];
ylabel('Share of total wealth (%)','fontsize',myfs)
xlabel('\tau_h','fontsize',myfs)
title(["(e) Top 1% wealth share",newline],'FontSize',myfs)

% legend

ax = nexttile(1);
lg = legend('Benchmark','Tax Reform','fontsize',myfs,'location','best');
lg.Layout.Tile = 6; % THIS IS THE COMMAND THAT WILL SET THE TILE WHERE YOU WANT TO PLOT YOUR LEGEND

if do_save==1; print([file_root,'SUERF'],'-dpng'); end
if do_save==1; print([file_root,'SUERF'],'-depsc'); end

% Up to here!
%% Make Plots: Figure 7
% Panel a: Total tax revenue (normalized)
% Panel b: Aggregate output (normalized)
% Panel c: Aggregate capital (norliazed)
% Panel d: Share of entre.
% Panel e: Entrepreneurial capital (normalized)
% Panel f: Share of entrepreneurial output
% Panel g: wage (normalized)
% Panel h: interest rate (normalized)
% 

% Panel a: Total tax revenue (normalized) 
[~,base_argmax] = max(base.taxes_inc_corp_div);
[~,CF2_argmax] = max(CF2.taxes_inc_corp_div);

%
%tau_h_grid = 1:numel(base.tau_h);
figure
plot(base.tau_h(tau_h_grid),base.taxes_inc_corp_div(tau_h_grid)/base.taxes_inc_corp_div(base_toptax_ind),...
    mylp1,'LineWidth',mylw+1,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),CF2.taxes_inc_corp_div(tau_h_grid)/CF2.taxes_inc_corp_div(base_toptax_ind),...
    mylp2,'LineWidth',mylw+1,'Color',mylc3)
hold on
xline(base.tau_h(base_argmax),mylp1,'LineWidth',mylw-1,'Color',mylc2)
hold on
xline(CF2.tau_h(CF2_argmax),mylp2,'LineWidth',mylw-1,'Color',mylc3)
%hold on
%xline(top_tax,'-k','LineWidth',mylw-1)
hold off
ax = gca;
ax.FontSize = myfs;
%ax.YLim = [0.85,1.2];
xlabel('\tau_h','fontsize',myfs)
ax.XLim=[base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))];
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
ylabel('1 = Baseline \tau_h','fontsize',myfs)
%title('Tax revenue ')
if do_save==1; print([file_root,'taxes_inc_corp_div'],'-dpng'); end
if do_save==1; print([file_root,'taxes_inc_corp_div'],'-depsc'); end

% Panel b: Aggregate output (normalized)
figure
plot(base.tau_h(tau_h_grid),(base.Y(tau_h_grid)/base.Y(base_toptax_ind)),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),(CF2.Y(tau_h_grid)/CF2.Y(CF2_toptax_ind)),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim = [0.85,1.2];
legend('Benchmark','Equal Tax Treatment','fontsize',myfs,'location','best')
xlabel('\tau_h','fontsize',myfs)
ylabel('1 = Baseline \tau_h')
if do_save==1; print([file_root,'Y'],'-dpng'); end
if do_save==1; print([file_root,'Y'],'-depsc'); end

% Panel c: Aggregate capital (norliazed)
figure
plot(base.tau_h(tau_h_grid),(base.K(tau_h_grid)/base.K(base_toptax_ind)),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),(CF2.K(tau_h_grid)/CF2.K(CF2_toptax_ind)),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim = [0.85,1.2];
% xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
ylabel('1 = Baseline \tau_h')
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'K'],'-dpng'); end
if do_save==1; print([file_root,'K'],'-depsc'); end

% Panel d: Share of entre.
figure
plot(base.tau_h(tau_h_grid),100*base.share_entre_act(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.share_entre_act(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim = [10, 20];
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
ylabel('Share of population (%)','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'share_entre_act'],'-dpng'); end
if do_save==1; print([file_root,'share_entre_act'],'-depsc'); end

% Panel e: Entrepreneurial capital (normalized)
figure
plot(base.tau_h(tau_h_grid),(base.aveK_all(tau_h_grid)/base.aveK_all(base_toptax_ind)),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),(CF2.aveK_all(tau_h_grid)/CF2.aveK_all(CF2_toptax_ind)),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim = [0.85,1.2];
% xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
ylabel('1 = Baseline \tau_h')
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'aveK_all'],'-dpng'); end
if do_save==1; print([file_root,'aveK_all'],'-depsc'); end

% Panel f: Share of entrepreneurial output
figure
plot(base.tau_h(tau_h_grid),100*base.Y_entre_share(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.Y_entre_share(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[66,76];
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
ylabel('Share of output (%)','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'Y_entre_share'],'-dpng'); end
if do_save==1; print([file_root,'Y_entre_share'],'-depsc'); end

% Panel g: wage (normalized)
figure
plot(base.tau_h(tau_h_grid),(base.w(tau_h_grid)/base.w(base_toptax_ind)),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),(CF2.w(tau_h_grid)/CF2.w(CF2_toptax_ind)),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim = [0.85,1.2];
% xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
ylabel('1 = Baseline \tau_h')
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'w'],'-dpng'); end
if do_save==1; print([file_root,'w'],'-depsc'); end


% Panel h: interest rate (normalized)
figure
plot(base.tau_h(tau_h_grid),(base.r(tau_h_grid)/base.r(base_toptax_ind)),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),(CF2.r(tau_h_grid)/CF2.r(CF2_toptax_ind)),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim = [0.85,1.2];
% xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
ylabel('1 = Baseline \tau_h')
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'r'],'-dpng'); end
if do_save==1; print([file_root,'r'],'-depsc'); end

%% Figure 8
% Panel a: Sole-prop share of entre
% Panel b: S-corp share of entre
% Panel c: C-corp share of entre
% Panel d: share of income benchmark
% Panel e: Taxes paid by workers
% Panel f: Taxes paid by sole prop
% Panel g: Taxes paid by S-corp
% Panel h: Taxes paid by C-corp

% Panel a: Sole-prop share of entre
figure
plot(base.tau_h(tau_h_grid),100*base.share_EP_entre(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.share_EP_entre(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[60,100];
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
ylabel('Share of entrepreneurs (%)','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%legend('Benchmark','Tax Reform','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'share_EP_entre'],'-dpng'); end
if do_save==1; print([file_root,'share_EP_entre'],'-depsc'); end


% Panel b: S-corp share of entre
figure
plot(base.tau_h(tau_h_grid),100*base.share_ES_entre(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.share_ES_entre(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
%hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[-1,30];
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
ylabel('Share of entrepreneurs (%)','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
legend('Benchmark','Tax Reform','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'share_ES_entre'],'-dpng'); end
if do_save==1; print([file_root,'share_ES_entre'],'-depsc'); end

% Panel c: C-corp share of entre
figure
plot(base.tau_h(tau_h_grid),100*base.share_EC_entre(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.share_EC_entre(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[-1,30];
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
ylabel('Share of entrepreneurs (%)','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'share_EC_entre'],'-dpng'); end
if do_save==1; print([file_root,'share_EC_entre'],'-depsc'); end


% Panel d: share of income declared as wage by S and C-corps, benchmark
figure
plot(base.tau_h(tau_h_grid),100*base.share_wage_ES(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(base.tau_h(tau_h_grid),100*base.share_wage_EC(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc2)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[5,45];
ylabel('Share of income (%)','fontsize',myfs)
%xline(top_tax,'-k','LineWidth',mylw-1)
legend('S-corp. (Benchmark)','C-corp. (Benchmark)','Location','best','fontsize',myfs)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%title('Share of income declared as wage')
if do_save==1; print([file_root,'base_share_wage'],'-dpng'); end
if do_save==1; print([file_root,'base_share_wage'],'-depsc'); end

% Panel e: Taxes (including soc sec) paid by workers
figure
plot(base.tau_h(tau_h_grid),100*base.taxes_occ(tau_h_grid,1)./base.taxes_young(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.taxes_occ(tau_h_grid,1)./CF2.taxes_young(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
xlabel('\tau_h','fontsize',myfs)
ax.XLim=[base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))];
ax.YLim=[60,65];
%legend('Benchmark','No Avoidance','fontsize',myfs-2,'location','best')
ylabel('Share of taxes revenue (%)','fontsize',myfs)
if do_save==1; print([file_root,'taxes_occ_w'],'-dpng'); end
if do_save==1; print([file_root,'taxes_occ_w'],'-depsc'); end

% Panel f: Taxes paid by sole prop
figure
plot(base.tau_h(tau_h_grid),100*base.taxes_occ(tau_h_grid,2)./base.taxes_young(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.taxes_occ(tau_h_grid,2)./CF2.taxes_young(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
ax.XLim=[base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))];
ax.YLim=[5,20];
%legend('Benchmark','No Avoidance','fontsize',myfs-2,'location','best')
ylabel('Share of taxes revenue (%)','fontsize',myfs)
if do_save==1; print([file_root,'taxes_occ_EP'],'-dpng'); end
if do_save==1; print([file_root,'taxes_occ_EP'],'-depsc'); end

% Panel g: Taxes paid by S-corp
figure
plot(base.tau_h(tau_h_grid),100*base.taxes_occ(tau_h_grid,3)./base.taxes_young(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.taxes_occ(tau_h_grid,3)./CF2.taxes_young(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
ax.XLim=[base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))];
ax.YLim=[-1,10];
%legend('Benchmark','No Avoidance','fontsize',myfs-2,'location','best')
ylabel('Share of taxes revenue (%)','fontsize',myfs)
if do_save==1; print([file_root,'taxes_occ_ES'],'-dpng'); end
if do_save==1; print([file_root,'taxes_occ_ES'],'-depsc'); end

% Panel h: Taxes paid by C-corp
figure
plot(base.tau_h(tau_h_grid),100*base.taxes_occ(tau_h_grid,4)./base.taxes_young(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.taxes_occ(tau_h_grid,4)./CF2.taxes_young(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
ax.XLim=[base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))];
ax.YLim=[15,30];
%legend('Benchmark','No Avoidance','fontsize',myfs-2,'location','best')
ylabel('Share of taxes revenue (%)','fontsize',myfs)
if do_save==1; print([file_root,'taxes_occ_EC'],'-dpng'); end
if do_save==1; print([file_root,'taxes_occ_EC'],'-depsc'); end

%% Figure 9
% Panel a: Income gini
% Panel b: Wealth gini
% Panel c: Top 1% income share
% Panel d: Top 1% Wealth share
% Panel e: Top 10% Income share
% Panel f: Top 10% Wealth share

% Panel a: Income gini
name = 'gini_inc_all';
figure
plot(base.tau_h(tau_h_grid),base.(name)(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),CF2.(name)(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ylabel('Gini coefficient','fontsize',myfs)
xlabel('\tau_h','fontsize',myfs)
ax.YLim=[0.54,0.58];
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%legend('Benchmark','Tax Reform','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,name],'-dpng'); end
if do_save==1; print([file_root,name],'-depsc'); end

% Panel b: Wealth gini
name = 'gini_wealth_all';
figure
plot(base.tau_h(tau_h_grid),base.(name)(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),CF2.(name)(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[0.80,0.84];
ylabel('Gini coefficient','fontsize',myfs)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
legend('Benchmark','Tax Reform','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,name],'-dpng'); end
if do_save==1; print([file_root,name],'-depsc'); end

% Panel c: Top 1% income share
name = 'inc_share_top1';
figure
plot(base.tau_h(tau_h_grid),100*base.(name)(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.(name)(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[12,18];
ylabel('Share of total income (%)','fontsize',myfs)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,name],'-dpng'); end
if do_save==1; print([file_root,name],'-depsc'); end

% Panel d: Top 1% wealth share
name = 'wealth_share_top1';
figure
plot(base.tau_h(tau_h_grid),100*base.(name)(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.(name)(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[17,23];
ylabel('Share of total wealth (%)','fontsize',myfs)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,name],'-dpng'); end
if do_save==1; print([file_root,name],'-depsc'); end


% Panel e: Top 10% income share
name = 'inc_share_top10';
figure
plot(base.tau_h(tau_h_grid),100*base.(name)(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.(name)(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[52,57];
ylabel('Share of total income (%)','fontsize',myfs)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,name],'-dpng'); end
if do_save==1; print([file_root,name],'-depsc'); end

% Panel f: Top 10% wealth share
name = 'wealth_share_top10';
figure
plot(base.tau_h(tau_h_grid),100*base.(name)(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),100*CF2.(name)(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc3)
hold on
ax = gca;
ax.FontSize = myfs;
ax.YLim=[63,68];
ylabel('Share of total wealth (%)','fontsize',myfs)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
%legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,name],'-dpng'); end
if do_save==1; print([file_root,name],'-depsc'); end



%% Other aggregate variables
%{

varNames = {'C','N','share_toptax','ss_cap','tau_p'};

for ii=1:numel(varNames)
    name = varNames{ii};
    figure
    plot(base.tau_h(tau_h_grid),(base.(name)(tau_h_grid)/base.(name)(base_toptax_ind)),'LineWidth',mylw)
    hold on
    plot(CF2.tau_h(tau_h_grid),(CF2.(name)(tau_h_grid)/CF2.(name)(CF2_toptax_ind)),'--','LineWidth',mylw)
    hold on
    ax = gca;
    ax.FontSize = myfs;
   % xline(top_tax,'-k','LineWidth',mylw-1)
    xlabel('\tau_h','fontsize',myfs)
    %xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
    ylabel('Normalized at baseline \tau_h')
    legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
    %title('Share of entrepreneurs')
    if do_save==1; print([file_root,name],'-dpng'); end
    if do_save==1; print([file_root,name],'-depsc'); end
end


%% Other plots
base.taxes_corp_div = base.taxes_corp + base.taxes_div;
figure
plot(base.tau_h(tau_h_grid),base.taxes_inc(tau_h_grid)/base.taxes_inc(base_toptax_ind),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(base.tau_h(tau_h_grid),base.taxes_corp_div(tau_h_grid)/base.taxes_corp_div(base_toptax_ind),...
    mylp2,'LineWidth',mylw,'Color',mylc2)
hold on
ax = gca;
ax.FontSize = myfs;
xlabel('\tau_h','fontsize',myfs)
xlim=([base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))]);
legend('Inc. tax','Corp. and div. tax','fontsize',myfs,'location','best')
ylabel('Normalized Tax Revenue','fontsize',myfs)
if do_save==1; print([file_root,'base_taxes'],'-dpng'); end
if do_save==1; print([file_root,'base_taxes'],'-depsc'); end

%% Share of taxes paid by each group

figure
plot(base.tau_h(tau_h_grid),base.taxes_occ(tau_h_grid,2)./base.taxes_young(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc2)
hold on
plot(base.tau_h(tau_h_grid),base.taxes_occ(tau_h_grid,3)./base.taxes_young(tau_h_grid),...
    mylp2,'LineWidth',mylw,'Color',mylc2)
hold on
plot(base.tau_h(tau_h_grid),base.taxes_occ(tau_h_grid,4)./base.taxes_young(tau_h_grid),...
    mylp3,'LineWidth',mylw,'Color',mylc2)
ax = gca;
ax.FontSize = myfs;
ax.XLim=[base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))];
ax.YLim=[0 0.3];
xlabel('\tau_h','fontsize',myfs)
legend('Sole-prop.','S-corp.','C-corp.','fontsize',myfs,'location','best')
ylabel('Share of Tax Revenue (Young Agents)','fontsize',myfs)
if do_save==1; print([file_root,'taxes_occ_benchmark'],'-dpng'); end
if do_save==1; print([file_root,'taxes_occ_benchmark'],'-depsc'); end

figure
plot(CF2.tau_h(tau_h_grid),CF2.taxes_occ(tau_h_grid,2)./CF2.taxes_young(tau_h_grid),...
    mylp1,'LineWidth',mylw,'Color',mylc3)
hold on
%plot(CF2.tau_h(tau_h_grid),CF2.taxes_occ(tau_h_grid,3)./CF2.taxes_young(tau_h_grid),...
%    mylp2,'LineWidth',mylw,'Color',mylc3)
%hold on
plot(CF2.tau_h(tau_h_grid),CF2.taxes_occ(tau_h_grid,4)./CF2.taxes_young(tau_h_grid),...
    mylp3,'LineWidth',mylw,'Color',mylc3)
hold off
ax = gca;
ax.FontSize = myfs;
ax.XLim=[base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))];
ax.YLim=[0 0.3];
xlabel('\tau_h','fontsize',myfs)
legend('Sole-prop.','C-corp.','fontsize',myfs-2,'location','best')
ylabel('Share of Tax Revenue (Young Agents)','fontsize',myfs)
if do_save==1; print([file_root,'taxes_occ_CF2'],'-dpng'); end
if do_save==1; print([file_root,'taxes_occ_CF2'],'-depsc'); end



%% Total taxes (including soc sec) paid by retired agents
figure
plot(base.tau_h(tau_h_grid),base.taxes_r(tau_h_grid)/base.taxes_r(1),'LineWidth',mylw)
hold on
plot(CF2.tau_h(tau_h_grid),CF2.taxes_r(tau_h_grid)/CF2.taxes_r(1),'--','LineWidth',mylw)
hold on
ax = gca;
ax.FontSize = myfs;
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
xlim=([base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))]);
legend('Benchmark','No Avoidance','fontsize',myfs-2,'location','best')
ylabel('Taxes paid by retirees','fontsize',myfs)
if do_save==1; print([file_root,'taxes_occ_r'],'-dpng'); end
if do_save==1; print([file_root,'taxes_occ_r'],'-depsc'); end


% CEV
if isfield(base,'cev')
    [~,base_argmax] = max(base.cev);
    [~,CF2_argmax] = max(CF2.cev);
    figure
    plot(base.tau_h,base.cev,'LineWidth',mylw+1)
    hold on
    plot(CF2.tau_h,CF2.cev,'LineWidth',mylw-1)
    hold on
    xline(base.tau_h(base_argmax),'--b','LineWidth',mylw+1)
    hold on
    xline(CF2.tau_h(CF2_argmax),'--r','LineWidth',mylw-1)
    hold on
    xline(top_tax,'-k','LineWidth',mylw-1)
    hold off
    ax = gca;
    ax.FontSize = myfs;
    legend('Benchmark','No Avoidance','Optimal (Benchmark)', ...
        'Optimal (No Avoidance)','Benchmark \tau_h','fontsize',myfs-1,'location','best')
    xlabel('\tau_h','fontsize',myfs)
    %xlim([base.tau_h(1) base.tau_h(end)])
    ylabel('CEV (%)','fontsize',myfs)
    if do_save==1; print([file_root,'cev_agg'],'-dpng'); end
    if do_save==1; print([file_root,'cev_agg'],'-depsc'); end

end

if isfield(base,'cev_aggcomp')
    [~,base_argmax] = max(base.cev);
    [~,CF2_argmax] = max(CF2.cev);
    figure
    plot(base.tau_h,base.cev_aggcomp,'LineWidth',mylw+1)
    hold on
    plot(CF2.tau_h,CF2.cev_aggcomp,'LineWidth',mylw-1)
    hold on
    xline(base.tau_h(base_argmax),'--b','LineWidth',mylw+1)
    hold on
    xline(CF2.tau_h(CF2_argmax),'--r','LineWidth',mylw-1)
    hold on
    xline(top_tax,'-k','LineWidth',mylw-1)
    hold off
    ax = gca;
    ax.FontSize = myfs;
    legend('Benchmark','No Avoidance','Optimal (Benchmark)', ...
        'Optimal (No Avoidance)','Benchmark \tau_h','fontsize',myfs-1,'location','best')
    xlabel('\tau_h','fontsize',myfs)
    %xlim([base.tau_h(1) base.tau_h(end)])
    ylabel('CEV (%)','fontsize',myfs)
    if do_save==1; print([file_root,'cev_aggcomp'],'-dpng'); end
    if do_save==1; print([file_root,'cev_aggcomp'],'-depsc'); end

end

if isfield(base,'cev_distcomp')
    [~,base_argmax] = max(base.cev);
    [~,CF2_argmax] = max(CF2.cev);
    figure
    plot(base.tau_h,base.cev_distcomp,'LineWidth',mylw+1)
    hold on
    plot(CF2.tau_h,CF2.cev_distcomp,'LineWidth',mylw-1)
    hold on
    xline(base.tau_h(base_argmax),'--b','LineWidth',mylw+1)
    hold on
    xline(CF2.tau_h(CF2_argmax),'--r','LineWidth',mylw-1)
    hold on
    xline(top_tax,'-k','LineWidth',mylw-1)
    hold off
    ax = gca;
    ax.FontSize = myfs;
    legend('Benchmark','No Avoidance','Optimal (Benchmark)', ...
        'Optimal (No Avoidance)','Benchmark \tau_h','fontsize',myfs-1,'location','best')
    xlabel('\tau_h','fontsize',myfs)
    %xlim([base.tau_h(1) base.tau_h(end)])
    ylabel('CEV (%)','fontsize',myfs)
    if do_save==1; print([file_root,'cev_distcomp'],'-dpng'); end
    if do_save==1; print([file_root,'cev_distcomp'],'-depsc'); end

end

% CEV by group (quartile and occupation)
if isfield(base,'cev_qo')
    figure
    b=bar(base.cev_qo,mybw);
    b(1).FaceColor = myfc1;
    b(2).FaceColor = myfc3;
    ax = gca;
    ax.FontSize = myfs;
    legend('Worker','Entrepreneur','fontsize',myfs,'location','best')
    ylim([barlb,barub])
    xlabel('Wealth quartiles','fontsize',myfs)
    ylabel('CEV (%)','fontsize',myfs)
    if do_save==1; print([file_root,'base_cev_qo'],'-dpng'); end
    if do_save==1; print([file_root,'base_cev_qo'],'-depsc'); end
end

if isfield(CF2,'cev_qo')
    figure
    b=bar(CF2.cev_qo,mybw);
    b(1).FaceColor = myfc1;
    b(2).FaceColor = myfc3;
    ax = gca;
    ax.FontSize = myfs;
    legend('Worker','Entrepreneur','fontsize',myfs,'location','best')
    ylim([barlb,barub])
    xlabel('Wealth quartiles','fontsize',myfs)
    ylabel('CEV (%)','fontsize',myfs)
    if do_save==1; print([file_root,'CF2_cev_qo'],'-dpng'); end
    if do_save==1; print([file_root,'CF2_cev_qo'],'-depsc'); end
end

% CEV by occupation and legal form
if isfield(base,'cev_z')
    figure
    X = categorical({'Worker','Sole-Prop.','S-Corp.','C-Corp.'});
    X = reordercats(X,{'Worker','Sole-Prop.','S-Corp.','C-Corp.'});
    b=bar(X',[base.cev_z,CF2.cev_z],mybw);
     b(1).FaceColor = myfc1;
    b(2).FaceColor = myfc3;
    ax = gca;
    ax.FontSize = myfs;
    legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
    ylim([barlb,barub])
    %xlabel('Wealth decile','fontsize',myfs)
    ylabel('CEV (%)','fontsize',myfs)
    if do_save==1; print([file_root,'cev_z'],'-dpng'); end
    if do_save==1; print([file_root,'cev_z'],'-depsc'); end
end

% CEV by wealth bins
if isfield(base,'cev_q')
    figure
    X = categorical({'Q1','Q2','Q3','Q4'});
    X = reordercats(X,{'Q1','Q2','Q3','Q4'});
    b=bar(X',[base.cev_q,CF2.cev_q],mybw);
    b(1).FaceColor = myfc1;
    b(2).FaceColor = myfc3;
    ax = gca;
    ax.FontSize = myfs;
    legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
    ylim([barlb,barub])
    xlabel('Wealth quartiles','fontsize',myfs)
    ylabel('CEV (%)','fontsize',myfs)
    if do_save==1; print([file_root,'cev_q'],'-dpng'); end
    if do_save==1; print([file_root,'cev_q'],'-depsc'); end
end


% Population shares



varNames = {'share_lfo_top1','share_lfo_top5','share_lfo_top10'};

for ii=1:numel(varNames)
    name = varNames{ii};
    figure
    plot(base.tau_h(tau_h_grid),base.(name)(tau_h_grid,1),'LineWidth',mylw)
    hold on
    plot(CF2.tau_h(tau_h_grid),CF2.(name)(tau_h_grid,1),'--','LineWidth',mylw)
    hold on
    ax = gca;
    ax.FontSize = myfs;
    %xline(top_tax,'-k','LineWidth',mylw-1)
    xlabel('\tau_h','fontsize',myfs)
    ylabel('Share of Sole Proprietors','fontsize',myfs)
    %xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
    legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
    %title('Share of entrepreneurs')
    if do_save==1; print([file_root,name,'_EP'],'-dpng'); end
    if do_save==1; print([file_root,name,'_EP'],'-depsc'); end

    figure
    plot(base.tau_h(tau_h_grid),base.(name)(tau_h_grid,2),'LineWidth',mylw)
    hold on
    plot(CF2.tau_h(tau_h_grid),CF2.(name)(tau_h_grid,2),'--','LineWidth',mylw)
    hold on
    ax = gca;
    ax.FontSize = myfs;
    %xline(top_tax,'-k','LineWidth',mylw-1)
    xlabel('\tau_h','fontsize',myfs)
    ylabel('Share of S-Corp.','fontsize',myfs)
    %xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
    legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
    %title('Share of entrepreneurs')
    if do_save==1; print([file_root,name,'_ES'],'-dpng'); end
    if do_save==1; print([file_root,name,'_ES'],'-depsc'); end

    figure
    plot(base.tau_h(tau_h_grid),base.(name)(tau_h_grid,3),'LineWidth',mylw)
    hold on
    plot(CF2.tau_h(tau_h_grid),CF2.(name)(tau_h_grid,3),'--','LineWidth',mylw)
    hold on
    ax = gca;
    ax.FontSize = myfs;
   % xline(top_tax,'-k','LineWidth',mylw-1)
    xlabel('\tau_h','fontsize',myfs)
    ylabel('Share of C-Corp.','fontsize',myfs)
    %xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
    legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
    %title('Share of entrepreneurs')
    if do_save==1; print([file_root,name,'_EC'],'-dpng'); end
    if do_save==1; print([file_root,name,'_EC'],'-depsc'); end

end







% Aggregate variables normlized at benchmark tau_h





% Interest rate in level
figure
plot(base.tau_h(tau_h_grid),base.r(tau_h_grid),'LineWidth',mylw)
hold on
plot(CF2.tau_h(tau_h_grid),CF2.r(tau_h_grid),'--','LineWidth',mylw)
hold on
ax = gca;
ax.FontSize = myfs;
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
ylabel('Interest rate','fontsize',myfs)
legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'r_level'],'-dpng'); end
if do_save==1; print([file_root,'r_level'],'-depsc'); end

% Share top tax bracket in level
figure
plot(base.tau_h(tau_h_grid),base.share_toptax(tau_h_grid),'LineWidth',mylw)
hold on
plot(CF2.tau_h(tau_h_grid),CF2.share_toptax(tau_h_grid),'--','LineWidth',mylw)
hold on
ax = gca;
ax.FontSize = myfs;
%xline(top_tax,'-k','LineWidth',mylw-1)
xlabel('\tau_h','fontsize',myfs)
%xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
ylabel('Share in top tax bracket','fontsize',myfs)
legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
%title('Share of entrepreneurs')
if do_save==1; print([file_root,'share_toptax_level'],'-dpng'); end
if do_save==1; print([file_root,'share_toptax_level'],'-depsc'); end

% Average capital and ability (theta) per entrepreneur
varNames = {'aveK_all','aveK_ep','aveK_es','aveK_ec',...
    'avetheta_all','avetheta_ep','avetheta_es','avetheta_ec'};

for ii=1:numel(varNames)
    name = varNames{ii};
    figure
    plot(base.tau_h(tau_h_grid),base.(name)(tau_h_grid)/base.(name)(base_toptax_ind),'LineWidth',mylw)
    hold on
    plot(CF2.tau_h(tau_h_grid),CF2.(name)(tau_h_grid)/CF2.(name)(CF2_toptax_ind),'--','LineWidth',mylw)
    hold on
    ax = gca;
    ax.FontSize = myfs;
   % xline(top_tax,'-k','LineWidth',mylw-1)
    xlabel('\tau_h','fontsize',myfs)
    %xlim([base.tau_h(1) base.tau_h(size(base.tau_h))])
    ylabel('Normalized at baseline \tau_h')
    legend('Benchmark','No Avoidance','fontsize',myfs,'location','best')
    %title('Share of entrepreneurs')
    if do_save==1; print([file_root,name],'-dpng'); end
    if do_save==1; print([file_root,name],'-depsc'); end

end

%}


end %end function "fun_plot_compstat"