%% Tax function
% Goal: plot average tax rate against multiples of average income

clear
clc
close all

global hsv0 hsv1 x_h tau_h

FS_pol = 14;

y_ave = 0.410725; % average income in the benchmark model
hsv0 = 0.855; % tax parameters in the benchmark model
hsv1 = 0.06;

x_h   = 6.5;
tau_h = 0.7;

xgrid = linspace(1,10,100)';
ygrid = zeros(length(xgrid),1);
netinc = zeros(length(xgrid),1);
for ii = 1:length(xgrid)
    ygrid(ii) = avetax(xgrid(ii));
    netinc(ii) = xgrid(ii)-xgrid(ii)*avetax(xgrid(ii));
end

figure
subplot(2,1,1)
plot(xgrid,ygrid,'LineWidth',3)
hold on
xlabel('Income','FontSize',FS_pol), axis tight
ylabel('Average income tax rate','FontSize',FS_pol)
grid on
legend('Benchmark: \tau = 0.06','Location','best','FontSize',FS_pol)
hold off

subplot(2,1,2)
plot(xgrid,xgrid,'--','LineWidth',3)
hold on
plot(xgrid,netinc,'LineWidth',3)
hold on
xlabel('Gross Income','FontSize',FS_pol), axis tight
ylabel('Net income','FontSize',FS_pol)
grid on
hold off

function [T] = tax(x)

global hsv0 hsv1 x_h tau_h

if x<x_h
    T = x-hsv0*(x).^(1-hsv1);
elseif x>=x_h
    T = x_h-hsv0*(x_h).^(1-hsv1)+tau_h*(x-x_h);
else
    error('smth wrong')
end


end

function F = avetax(x)

F = tax(x)./x;

end
