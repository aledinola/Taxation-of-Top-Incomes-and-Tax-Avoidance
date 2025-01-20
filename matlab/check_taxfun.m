%% Check income tax function with top rate for rich
% Points where VFI does not converge:
% % 9
%     12
%     13
%     35
%     47
clear
clc
close all


x_min = 0.392;
x_max = 0.52;
n_points = 31;
x_int = (x_max-x_min)/(n_points-1);

n = round((0.6-x_min)/x_int+1)

% top_rate = 0.369;
% n_tau_h = 31;
% linspace(top_rate,0.53,n_tau_h)';
% 
% 
% n_y = 100;
% tau_H_values = [0.396,0.42]';
% mtr_vec      = zeros(n_y,length(tau_H_values));
% 
% for i_H=1:length(tau_H_values)
% 
% tau    = 0.127;
% tau_H  = tau_H_values(i_H);
% lambda = 0.796;
% y_H    = 3.01614179329491;%fun_y_H(lambda,tau,tau_H);
% 
% y_vec = linspace(0.05,1.5*y_H,n_y)';
% 
% mtr_vec(:,i_H) = fun_mtr(y_vec,lambda,tau,y_H,tau_H);
% 
% end
% 
% figure
% plot(y_vec,mtr_vec(:,1),'-o')
% hold on
% plot(y_vec,mtr_vec(:,2))
% xline(y_H,'--','LineWidth',2)
% legend('Bench','Higher tau_H')
% ylim([-inf,0.6])
% title('Marginal income tax rate')
% xlabel('Income')
% ylabel('Marginal rate')
% 
% function T = fun_tax(y,lambda,tau,y_H,tau_H)
% 
%     if y<=y_H
%         T = y-lambda*y.^(1-tau);
%     else
%         T = y_H-lambda*y_H^(1-tau)+tau_H*(y-y_H);
%     end
% 
% end
% 
% function T = fun_mtr(y,lambda,tau,y_H,tau_H)
% % Marginal tax rate
% 
%     n = length(y );
%     T = zeros(n,1);
%     for ii=1:n
%         y_val = y(ii);
%         if y_val<=y_H
%             T(ii) = 1-lambda*(1-tau)*y_val^(-tau);
%         else
%             T(ii) = tau_H;
%         end
%     end
% 
% end
% 
% function y_H = fun_y_H(lambda,tau,tau_H)
%     y_H = (lambda*(1-tau)/(1-tau_H))^(1/tau);
% end
% 
