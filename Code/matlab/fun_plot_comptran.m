function fun_plot_comptran(base,CF2,ss_base,ss_CF2, file_root,do_save)
% PURPOSE
% fun_plot_comptran makes plots of results base,CF2,ss_base,ss_CF2 in files 
% with stem filestem
% INPUTS
% base,CF2  : struct containing results from welfare maximization exercise
%             with transition, under benchmark and no-avoidance (cf2).
% ss_base,ss_CF2: struct containing initial steady state results, for
%                 benchmark and no-avoidance (cf2)
% file_root  : character with root of the filename
% do_save   : flag 0/1, if you want to save the plots

tvec = (1:size(base.KN,1))'; % Vector with time periods in transition
tau_h_bench = 0.396;

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


tau_h_grid = 1:38;
% CEV agg
[~,base_argmax] = max(base.cev_vec);
[~,CF2_argmax] = max(CF2.cev_vec);

disp(base.tau_h(tau_h_grid(base_argmax)))
disp(CF2.tau_h(tau_h_grid(CF2_argmax)))


% Aggregate CEV

% convert to percent
base.cev_vec = 100*base.cev_vec;
CF2.cev_vec = 100*CF2.cev_vec;
base.cev_qo = 100*base.cev_qo;
CF2.cev_qo = 100*CF2.cev_qo;
base.cev_z = 100*base.cev_z;
CF2.cev_z = 100*CF2.cev_z;

base.cev_aggcomp_vec = 100*base.cev_aggcomp_vec;
CF2.cev_aggcomp_vec = 100*CF2.cev_aggcomp_vec;
base.cev_distcomp_vec = 100*base.cev_distcomp_vec;
CF2.cev_distcomp_vec = 100*CF2.cev_distcomp_vec;



%% Make plots: Figure 11
% Panel a: Aggregate CEV
% Panel b: Agg component of CEV
% Panel c: Dist component of CEV
% Panel d: CEV by occ and LFO
% Panel e: CEV by wealth and occ, benchmark
% Panel f: CEV by wealth and occ, CF avoidance

% Panel a: Aggregate CEV
figure
plot(base.tau_h(tau_h_grid),base.cev_vec(tau_h_grid),...
    mylp1,'LineWidth',mylw+1,'Color',mylc2)
hold on
plot(CF2.tau_h(tau_h_grid),CF2.cev_vec(tau_h_grid),...
    mylp2,'LineWidth',mylw+1,'Color',mylc3)
hold on
xline(base.tau_h(base_argmax),mylp1,'LineWidth',mylw-1,'Color',mylc2)
hold on
xline(CF2.tau_h(CF2_argmax),mylp2,'LineWidth',mylw-1,'Color',mylc3)
ax = gca;
ax.FontSize = myfs;
xlabel('\tau_h','fontsize',myfs)
ax.XLim=([base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))]);
ylabel('CEV (%)','fontsize',myfs)
legend('Benchmark','Equal Tax Treatment','fontsize',myfs-1,'location','best')
if do_save==1; print([file_root,'cev_vec'],'-dpng'); end
if do_save==1; print([file_root,'cev_vec'],'-depsc'); end

% Decomposition of CEV
if isfield(base,'cev_aggcomp_vec')
    % Panel b: Agg component of CEV
    figure
    plot(base.tau_h(tau_h_grid),base.cev_aggcomp_vec(tau_h_grid)-base.cev_aggcomp_vec(1),...
        mylp1,'LineWidth',mylw,'Color',mylc2)
    hold on
    plot(CF2.tau_h(tau_h_grid),CF2.cev_aggcomp_vec(tau_h_grid)-CF2.cev_aggcomp_vec(1),...
        mylp2,'LineWidth',mylw,'Color',mylc3)   
    ax = gca;
    ax.FontSize = myfs;
    xlabel('\tau_h','fontsize',myfs)
    xlim=([base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))]);
    ylabel('Agg. component of CEV (%)','fontsize',myfs)
    %legend('Benchmark','No Avoidance','fontsize',myfs-1,'location','best')
    if do_save==1; print([file_root,'cev_aggcomp_vec'],'-dpng'); end
    if do_save==1; print([file_root,'cev_aggcomp_vec'],'-depsc'); end

% Panel c: Dist component of CEV
    figure
    plot(base.tau_h(tau_h_grid),base.cev_distcomp_vec(tau_h_grid)-base.cev_distcomp_vec(1),...
        mylp1,'LineWidth',mylw,'Color',mylc2)
    hold on
    plot(CF2.tau_h(tau_h_grid),CF2.cev_distcomp_vec(tau_h_grid)-CF2.cev_distcomp_vec(1),...
        mylp2,'LineWidth',mylw,'Color',mylc3)
    ax = gca;
    ax.FontSize = myfs;
    xlabel('\tau_h','fontsize',myfs)
    xlim=([base.tau_h(tau_h_grid(1)),base.tau_h(tau_h_grid(numel(tau_h_grid)))]);
    ylabel('Dist. component of CEV (%)','fontsize',myfs)
    %legend('Benchmark','No Avoidance','fontsize',myfs-1,'location','best')
    if do_save==1; print([file_root,'cev_distcomp_vec'],'-dpng'); end
    if do_save==1; print([file_root,'cev_distcomp_vec'],'-depsc'); end
end

% Panel d: CEV by occ and LFO
CF2.cev_z(3)=0;
figure
X = categorical({'Worker','Sole-Prop.','S-Corp.','C-Corp.'});
X = reordercats(X,{'Worker','Sole-Prop.','S-Corp.','C-Corp.'});
b=bar(X',[base.cev_z,CF2.cev_z],mybw,'LineWidth',mybaredge);
b(1).FaceColor = myfc1;
b(2).FaceColor = myfc3;
%b(3).FaceColor = myfc3;
ax = gca;
ax.FontSize = myfs;
legend('Benchmark','Equal Tax Treatment','fontsize',myfs,'location','southwest')
%ylim([-6e-3,1e-2])
ylabel('CEV (%)','fontsize',myfs)
if do_save==1; print([file_root,'cev_z'],'-dpng'); end
if do_save==1; print([file_root,'cev_z'],'-depsc'); end

% Panel e: CEV by wealth and occ, benchmark
figure
b=bar(base.cev_qo,mybw,'LineWidth',mybaredge);
b(1).FaceColor = myfc1;
b(2).FaceColor = myfc3;
ax = gca;
ax.FontSize = myfs;
legend('Worker','Entrepreneur','fontsize',myfs,'location','best')
ylim([-6e-1,1e0])
xlabel('Wealth quartile','fontsize',myfs)
ylabel('CEV (%)','fontsize',myfs)
if do_save==1; print([file_root,'base_cev_qo'],'-dpng'); end
if do_save==1; print([file_root,'base_cev_qo'],'-depsc'); end

% Panel f: CEV by wealth and occ, CF avoidance
figure
b=bar(CF2.cev_qo,mybw,'LineWidth',mybaredge);
b(1).FaceColor = myfc1;
b(2).FaceColor = myfc3;
legend('Worker','Entrepreneur','fontsize',myfs,'location','best')
ax = gca;
ax.FontSize = myfs;
ylim([-6e-1,1e0])
xlabel('Wealth quartile','fontsize',myfs)
ylabel('CEV (%)','fontsize',myfs)
if do_save==1; print([file_root,'CF2_cev_qo'],'-dpng'); end
if do_save==1; print([file_root,'CF2_cev_qo'],'-depsc'); end

end %end function fun_plot_comptran