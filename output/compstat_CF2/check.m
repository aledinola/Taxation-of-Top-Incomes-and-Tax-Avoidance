% compstat_CF2
clear
clc
close all

% "taxes_occ" is a (n_tau_h,no) array, with no=4, which gives the total taxes
% paid by each of the no groups of young people: {W,EP,ES,EC}. 
% Note that in CF2 (no avoidance), the ES group is empty, so
% taxes_occ(tau_h,3) = 0, for all tau_h.
% "taxes_tot" gives total taxes, including social security taxes

tau_h_vec  = load('tau_h_vec.txt');

taxes_tot  = load('taxes_tot_vec.txt');
taxes_corp = load('taxes_corp_vec.txt');
taxes_div  = load('taxes_div_vec.txt');
taxes_inc  = load('taxes_inc_vec.txt');
taxes_ss   = load('taxes_ss_vec.txt');
taxes_occ  = load('taxes_occ_vec.txt');
taxes_r    = load('taxes_r_vec.txt');
Y_vec      = load('Y_vec.txt');
share_EC_entre = load('share_EC_entre_vec.txt');
taxes_occ  = reshape(taxes_occ,[numel(tau_h_vec),4]);
taxes_inc_corp_div = taxes_inc + taxes_div + taxes_corp;

[~,max_ind] = max(taxes_inc_corp_div);
tau_h_max = tau_h_vec(max_ind);
disp('Revenue-maximizing tax rate')
disp(tau_h_max)

% Check
check_taxes = sum(taxes_occ,2)+taxes_r-taxes_tot;
if any(abs(check_taxes)>1e-6)
    warning('taxes tot is wrong')
end

% Make plots
tau_h_grid = 1:numel(tau_h_vec);

figure
plot(tau_h_vec(tau_h_grid),taxes_inc_corp_div(tau_h_grid)/taxes_inc_corp_div(1))
title('taxes inc + corp + div')

figure
plot(tau_h_vec(tau_h_grid),taxes_tot(tau_h_grid)/taxes_tot(1))
title('taxes inc + corp + div + socsec')

figure
plot(tau_h_vec,taxes_ss/taxes_ss(1))
title('Social security taxes')

figure
plot(tau_h_vec,taxes_occ(:,1))
title('taxes paid by W')
figure
plot(tau_h_vec,taxes_occ(:,2))
title('taxes paid by EP')
figure
plot(tau_h_vec,taxes_occ(:,3))
title('taxes paid by ES')
figure
plot(tau_h_vec,taxes_occ(:,4))
title('taxes paid by EC')

% n = 32;
% 
% x_min = 0.396;
% x_max = 0.52;
% step = (x_max-x_min)/(n-1);
% 
% % Now we set x_max = 0.55
% % Keep step constant, what is the right n?
% x_max = 0.552;
% n1 = (x_max-x_min)/step+1;