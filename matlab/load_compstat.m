function A = load_compstat(dirname)
% PURPOSE
%load_compstat loads comparative statics result from file 'filename' and
%stores them into the structure A
% INPUTS
% filename : name of the txt directory where the results of comparative statics 
%            are saved
% OUTPUTS
% A        : structure where variables are stored

if ~ischar(dirname)
    error('Input argument dirname must be a character')
end

%result = load(dirname);

% Read the size of the first dimension
tmp = load(fullfile(dirname,'tau_h_vec.txt'));
n_tau_h = numel(tmp);

siz1 = [n_tau_h,1];
siz3 = [n_tau_h,3];
siz4 = [n_tau_h,4];

nq = 4;nz = 4;

%% Fields are 1-dim, with size = [n_tau_h,1]
A.tau_h        = loadArray(fullfile(dirname,'tau_h_vec.txt'),siz1);
A.tau_p        = loadArray(fullfile(dirname,'tau_p0_vec.txt'),siz1);
A.hsv_0        = loadArray(fullfile(dirname,'hsv_0_vec.txt'),siz1);
A.taxes_tot    = loadArray(fullfile(dirname,'taxes_tot_vec.txt'),siz1);
A.taxes_ss     = loadArray(fullfile(dirname,'taxes_ss_vec.txt'),siz1);
A.taxes_inc    = loadArray(fullfile(dirname,'taxes_inc_vec.txt'),siz1);
A.taxes_corp   = loadArray(fullfile(dirname,'taxes_corp_vec.txt'),siz1);
A.taxes_div    = loadArray(fullfile(dirname,'taxes_div_vec.txt'),siz1);
A.res_gov      = loadArray(fullfile(dirname,'res_gov_vec.txt'),siz1);
A.share_toptax = loadArray(fullfile(dirname,'share_toptax_vec.txt'),siz1);
A.Y = loadArray(fullfile(dirname,'Y_vec.txt'),siz1);
A.Y_entre = loadArray(fullfile(dirname,'Y_entre_vec.txt'),siz1);
A.Y_entre_share = A.Y_entre/A.Y;
A.r = loadArray(fullfile(dirname,'r_vec.txt'),siz1);
A.w = loadArray(fullfile(dirname,'w_vec.txt'),siz1);
A.ss_cap = loadArray(fullfile(dirname,'ss_cap_vec.txt'),siz1);

A.share_entre_act   = loadArray(fullfile(dirname,'share_entre_act_vec.txt'),siz1);
A.share_EP_entre    = loadArray(fullfile(dirname,'share_EP_entre_vec.txt'),siz1);
A.share_ES_entre    = loadArray(fullfile(dirname,'share_ES_entre_vec.txt'),siz1);
A.share_EC_entre    = loadArray(fullfile(dirname,'share_EC_entre_vec.txt'),siz1);
A.gini_wealth_all   = loadArray(fullfile(dirname,'gini_wealth_all_vec.txt'),siz1);
A.gini_wealth_entre = loadArray(fullfile(dirname,'gini_wealth_entre_vec.txt'),siz1);
A.gini_inc_all      = loadArray(fullfile(dirname,'gini_inc_all_vec.txt'),siz1);
A.gini_inc_entre    = loadArray(fullfile(dirname,'gini_inc_entre_vec.txt'),siz1);
A.K_entre           =  loadArray(fullfile(dirname,'K_entre_vec.txt'),siz1);
A.K_C               =  loadArray(fullfile(dirname,'K_C_vec.txt'),siz1);
A.K                 =  A.K_C + A.K_entre;
A.N_entre           =  loadArray(fullfile(dirname,'N_entre_vec.txt'),siz1);
A.N_C               =  loadArray(fullfile(dirname,'N_C_vec.txt'),siz1);
A.N                 =  A.N_C + A.N_entre;
A.C                 =  loadArray(fullfile(dirname,'C_vec.txt'),siz1);
A.share_wage_ES     =  loadArray(fullfile(dirname,'share_wage_ES_vec.txt'),siz1);
A.share_wage_EC     =  loadArray(fullfile(dirname,'share_wage_EC_vec.txt'),siz1);

%% Fields are 2-dim, with size = [n_tau_h,no], where no=4
A.wealth_share          = loadArray(fullfile(dirname,'wealth_share_vec.txt'),siz4);
A.wealth_share_top1     = A.wealth_share(:,1);
A.wealth_share_top10    = A.wealth_share(:,2);
A.wealth_share_top20    = A.wealth_share(:,3);
A.wealth_share_bottom40 = A.wealth_share(:,4);

A.inc_share          = loadArray(fullfile(dirname,'inc_share_vec.txt'),siz4);
A.inc_share_top1     = A.inc_share(:,1);
A.inc_share_top10    = A.inc_share(:,2);
A.inc_share_top20    = A.inc_share(:,3);
A.inc_share_bottom40 = A.inc_share(:,4);

A.aveK     = loadArray(fullfile(dirname,'avek_vec.txt'),siz4);
A.aveK_all = A.aveK(:,1); %all entre
A.aveK_ep  = A.aveK(:,2); %sole-prop 
A.aveK_es  = A.aveK(:,3); %S-corp
A.aveK_ec  = A.aveK(:,4); %C-corp

A.avetheta = loadArray(fullfile(dirname,'avetheta_vec.txt'),siz4);
A.avetheta_all = A.avetheta(:,1); %all entre
A.avetheta_ep  = A.avetheta(:,2); %sole-prop 
A.avetheta_es  = A.avetheta(:,3); %S-corp
A.avetheta_ec  = A.avetheta(:,4); %C-corp


A.taxes_top = loadArray(fullfile(dirname,'taxes_top_vec.txt'),siz4);
A.taxes_top1  = A.taxes_top(:,1);
A.taxes_top3  = A.taxes_top(:,2);
A.taxes_top5  = A.taxes_top(:,3);
A.taxes_top10 = A.taxes_top(:,4);

% NOTE
% "taxes_occ" is a (n_tau_h,no) array, with no=4, which gives the total taxes
% paid by each of the no groups of young people: {W,EP,ES,EC}. 
% Note that in CF2 (no avoidance), the ES group is empty, so
% taxes_occ(tau_h,3) = 0, for all tau_h.
% "taxes_tot" gives total taxes, including social security taxes
A.taxes_occ = loadArray(fullfile(dirname,'taxes_occ_vec.txt'),siz4);
A.taxes_r   = loadArray(fullfile(dirname,'taxes_r_vec.txt'),siz1);

%% Fields are 2-dim, with size = [n_tau_h,3]
A.share_lfo_top1  = loadArray(fullfile(dirname,'share_lfo_top1.txt'),siz3);
A.share_lfo_top5  = loadArray(fullfile(dirname,'share_lfo_top5.txt'),siz3);
A.share_lfo_top10 = loadArray(fullfile(dirname,'share_lfo_top10.txt'),siz3);


%% CEV (only for fiscal neutrality case)
if isfile(fullfile(dirname,'cev_vec.txt'))
    A.cev = 100*loadArray(fullfile(dirname,'cev_vec.txt'),siz1);
end
if isfile(fullfile(dirname,'cev_aggcomp_vec.txt'))
    A.cev_aggcomp = 100*loadArray(fullfile(dirname,'cev_aggcomp_vec.txt'),siz1);
end
if isfile(fullfile(dirname,'cev_distcomp_vec.txt'))
    A.cev_distcomp = 100*loadArray(fullfile(dirname,'cev_distcomp_vec.txt'),siz1);
end
if isfile(fullfile(dirname,'cev_q.txt'))
    A.cev_q = 100*loadArray(fullfile(dirname,'cev_q.txt'),[nq,1]);
end
if isfile(fullfile(dirname,'cev_qo.txt'))
    A.cev_qo = 100*loadArray(fullfile(dirname,'cev_qo.txt'),[nq,2]);
end
if isfile(fullfile(dirname,'cev_z.txt'))
    A.cev_z = 100*loadArray(fullfile(dirname,'cev_z.txt'),[nz,1]);
end

%% Create additional variables
A.taxes_corp_div     = A.taxes_corp + A.taxes_div;
A.taxes_inc_corp_div = A.taxes_inc + A.taxes_corp + A.taxes_div;
A.Y_net              = A.Y - A.taxes_tot;
%A.taxes_tot_tr       = A.taxes_tot - A.tr;

if ~isstruct(A)
    error('Output argument A must be a structure')
end

end %end function "load_compstat"