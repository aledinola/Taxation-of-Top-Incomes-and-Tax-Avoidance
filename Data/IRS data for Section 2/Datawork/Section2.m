
clear;
clc;

load data;       % IRS data from SOI
load wageshare;  % selected data by business receipts, S-corp and C-corp for the years 2007, 2013, 2019

% Statistics from Integrated Business Data

time_IBD = 1980:1:2015;
le = size(time_IBD);
maxle = le(1,2);


share_corp = 100*IRSIBD(7,:)./IRSIBD(1,:); % corporation/total businesses
share_Ccorp = 100*IRSIBD(13,:)./IRSIBD(1,:); % C-corporation/total businesses
share_Scorp = 100*IRSIBD(25,:)./IRSIBD(1,:); % S-corporation/total businesses
share_partner = 100*IRSIBD(31,:)./IRSIBD(1,:); % partnerships/total businesses
share_sole = 100*IRSIBD(55,:)./IRSIBD(1,:);  % sole props/total businesses

share_passthroughs = share_Scorp + share_sole + share_partner;

fig1 = plot(time_IBD, share_passthroughs,'k-', time_IBD, share_Scorp, 'b--', time_IBD, share_Ccorp, 'r-.')
xlabel(' Year');
ylabel(' Share of businesses in %');
%title('\bf Share of Businesses')
%legend('Passthroughs', 'S-corporations', 'C-corporations', 'Location', 'NorthWest');
axis([time_IBD(1) time_IBD(maxle) -10 110]);
set(fig1,  'LineWidth',3);
set(gca, 'FontSize',18);
%saveas(gca, 'shareLFOversion1','epsc')

%pause;
close all;

yyaxis right
plot(time_IBD, share_passthroughs, 'k-', 'LineWidth', 3);
ylim([70 100])
ylabel('Share of pass-throughs in %', 'Rotation', 270);
set(gca, 'FontSize', 18, 'YColor', 'k');
yyaxis left
plot(time_IBD, share_Scorp, 'b--', time_IBD, share_Ccorp, 'r-.', 'LineWidth', 3);
ylim([0 20])
ylabel(' Share of S-corp., C-corp. in %');
xlim([time_IBD(1), time_IBD(maxle)])
xlabel('Year');
%title('Share of Businesses');
set(gca, 'FontSize', 18, 'YColor', 'k');
saveas(gca, 'shareLFO', 'epsc');

%pause;
close all;


share_corp_netincome = 100*IRSIBD(10,:)./IRSIBD(4,:);  % net income less deficit of corporations / net income less deficit of all businesses
share_Ccorp_netincome = 100*IRSIBD(16,:)./IRSIBD(4,:);  % net income less deficit of C-corporations / net income less deficit of all businesses
share_Scorp_netincome = 100*IRSIBD(28,:)./IRSIBD(4,:);  % net income less deficit of S-corporations / net income less deficit of all businesses
share_partner_netincome = 100*IRSIBD(34,:)./IRSIBD(4,:);  % net income less deficit of partnerships / net income less deficit of all businesses
share_sole_netincome = 100*IRSIBD(58,:)./IRSIBD(4,:);  % net income less deficit of sole proprietors / net income less deficit of all businesses
share_passthroughs_netincome = share_Scorp_netincome + share_sole_netincome + share_partner_netincome;


fig2 = plot(time_IBD, share_passthroughs_netincome,'k-', time_IBD, share_Scorp_netincome, 'b--', time_IBD, share_Ccorp_netincome, 'r-.')
xlabel('Year');
ylabel(' Share of net income in %');
%title('\bf Share of Net Income (less deficit)')
legend('Pass-throughs', 'S-corporations', 'C-corporations', 'Location', 'NorthEast');
axis([time_IBD(1) time_IBD(maxle) -10 110]);
set(fig2,  'LineWidth',3);
set(gca, 'FontSize',18);
saveas(gca, 'incomeLFO','epsc')

%pause;
close all;

share_corp_busreceipts = 100*IRSIBD(9,:)./IRSIBD(3,:);  % net income less deficit of corporations / net income less deficit of all businesses
share_Ccorp_busreceipts = 100*IRSIBD(15,:)./IRSIBD(3,:);  % net income less deficit of C-corporations / net income less deficit of all businesses
share_Scorp_busreceipts = 100*IRSIBD(27,:)./IRSIBD(3,:);  % net income less deficit of S-corporations / net income less deficit of all businesses
share_partner_busreceipts = 100*IRSIBD(33,:)./IRSIBD(3,:);  % net income less deficit of partnerships / net income less deficit of all businesses
share_sole_busreceipts = 100*IRSIBD(57,:)./IRSIBD(3,:);  % net income less deficit of sole proprietors / net income less deficit of all businesses
share_passthroughs_busreceipts = share_Scorp_busreceipts + share_sole_busreceipts + share_partner_busreceipts;


fig2 = plot(time_IBD, share_passthroughs_busreceipts,'k-', time_IBD, share_Scorp_busreceipts, 'b--', time_IBD, share_Ccorp_busreceipts, 'r-.')
xlabel('Year');
ylabel(' Share of business receipts in %');
%title('\bf Share of Net Income (less deficit)')
legend('Pass-throughs', 'S-corporations', 'C-corporations', 'Location', 'NorthEast');
axis([time_IBD(1) time_IBD(maxle) -10 110]);
set(fig2,  'LineWidth',3);
set(gca, 'FontSize',18);
saveas(gca, 'busreceiptsLFO','epsc')


close all;



% Statistics from IRS

time_IRS = 1995:1:2020;
le = size(time_IRS);
maxle = le(1,2);

C_compensation = CCorp(47,:);  % officer compensation, C-corporation
C_netincome    = CCorp(65,:);  % net income less deficit, C-corporation
C_wage_share = 100.*C_compensation./(C_compensation + C_netincome);


S_compensation = SCorp(10,:);  % officer compensation, C-corporation
S_netincome    = SCorp(42,:);  % net income less deficit, C-corporation
S_wage_share = 100.*S_compensation./(S_compensation + S_netincome);


% fig3 = plot(time_IRS, S_wage_share,'b--', time_IRS, C_wage_share, 'r-.')
% xlabel('year');
% ylabel(' %');
% %title('\bf Business Owner Wages as Share of Net Income')
% %egend('S-corporations', 'C-corporations', 'Location', 'NorthWest');
% axis([time_IRS(1) time_IRS(maxle) 5 50]);
% set(fig3,  'LineWidth',3);
% set(gca, 'FontSize',18);
% saveas(gca, 'wageLFO','epsc')
% 
% %pause;
% close all;

fig3 = plot(time_IRS, S_wage_share,'b--')
xlabel('Year');
ylabel(' Wage share in % ');
%title('\bf Business Owner Wages as Share of Net Income')
%egend('S-corporations', 'C-corporations', 'Location', 'NorthWest');
axis([time_IRS(1) time_IRS(maxle) 20 60]);
set(fig3,  'LineWidth',3);
set(gca, 'FontSize',18);
saveas(gca, 'wageLFO','epsc')

%pause;
close all;


% statistics from NIPA

time_NIPA = 1971:1:2023;
le = size(time_NIPA);
maxle = le(1,2);

laborshare = 100.*NIPA(4,:)./NIPA(1,:); % compensation of employees/gross value added of corporate businesses.

fig3b = plot(time_NIPA, laborshare,'b-')
xlabel(' Year');
ylabel(' Corporate labor share in %');
%title('\bf Labor Share in Corporate Sector')
axis([time_NIPA(1) time_NIPA(maxle) 50 70]);
set(fig3b,  'LineWidth',3);
set(gca, 'FontSize',18);
saveas(gca, 'laborshare','epsc')


%pause;
close all;


% statistics on wage share for S-corporations and C-corporations

wage_Scorp_2007 = Swage2007(4,:)./(Swage2007(4,:)+Swage2007(6,:));
wage_Scorp_2013 = Swage2013(4,:)./(Swage2013(4,:)+Swage2013(6,:));
wage_Scorp_2019 = Swage2019(4,:)./(Swage2019(4,:)+Swage2019(6,:));

wage_Ccorp_2007 = Cwage2007(4,:)./(Cwage2007(4,:)+Cwage2007(6,:));
wage_Ccorp_2013 = Cwage2013(4,:)./(Cwage2013(4,:)+Cwage2013(6,:));
wage_Ccorp_2019 = Cwage2019(4,:)./(Cwage2019(4,:)+Cwage2019(6,:));


histo_Scorp_2007 = [100.*sum(Swage2007(4,2:6),2)./(sum(Swage2007(4,2:6),2)+sum(Swage2007(6,2:6),2)) % business receipts less than 1. Mio
                    100.*sum(Swage2007(4,7:8),2)./(sum(Swage2007(4,7:8),2)+sum(Swage2007(6,7:8),2)) % business receipt between 1 and 5 Mio
                    100.*sum(Swage2007(4,9),2)./(sum(Swage2007(4,9),2)+sum(Swage2007(6,9),2)) % business receipt between 5 and 10 Mio
                    100.*sum(Swage2007(4,10:11),2)./(sum(Swage2007(4,10:11),2)+sum(Swage2007(6,10:11),2))]; % business receipt larger than 10 Mio

histo_Ccorp_2007 = [100.*sum(Cwage2007(4,2:6),2)./(sum(Cwage2007(4,2:6),2)+sum(Cwage2007(6,2:6),2)) % business receipts less than 1. Mio
                    100.*sum(Cwage2007(4,7:8),2)./(sum(Cwage2007(4,7:8),2)+sum(Cwage2007(6,7:8),2)) % business receipt between 1 and 5 Mio
                    100.*sum(Cwage2007(4,9),2)./(sum(Cwage2007(4,9),2)+sum(Cwage2007(6,9),2)) % business receipt between 5 and 10 Mio
                    100.*sum(Cwage2007(4,10:13),2)./(sum(Cwage2007(4,10:13),2)+sum(Cwage2007(6,10:13),2))]; % business receipt larger than 10 Mio

histo_Scorp_2013 = [100.*sum(Swage2013(4,2:6),2)./(sum(Swage2013(4,2:6),2)+sum(Swage2013(6,2:6),2)) % business receipts less than 1. Mio
                    100.*sum(Swage2013(4,7:8),2)./(sum(Swage2013(4,7:8),2)+sum(Swage2013(6,7:8),2)) % business receipt between 1 and 5 Mio
                    100.*sum(Swage2013(4,9),2)./(sum(Swage2013(4,9),2)+sum(Swage2013(6,9),2)) % business receipt between 5 and 10 Mio
                    100.*sum(Swage2013(4,10:11),2)./(sum(Swage2013(4,10:11),2)+sum(Swage2013(6,10:11),2))]; % business receipt larger than 10 Mio

histo_Ccorp_2013 = [100.*sum(Cwage2013(4,2:6),2)./(sum(Cwage2013(4,2:6),2)+sum(Cwage2013(6,2:6),2)) % business receipts less than 1. Mio
                    100.*sum(Cwage2013(4,7:8),2)./(sum(Cwage2013(4,7:8),2)+sum(Cwage2013(6,7:8),2)) % business receipt between 1 and 5 Mio
                    100.*sum(Cwage2013(4,9),2)./(sum(Cwage2013(4,9),2)+sum(Cwage2013(6,9),2)) % business receipt between 5 and 10 Mio
                    100.*sum(Cwage2013(4,10:13),2)./(sum(Cwage2013(4,10:13),2)+sum(Cwage2013(6,10:13),2))]; % business receipt larger than 10 Mio



histo_Scorp_2019 = [100.*sum(Swage2019(4,2:6),2)./(sum(Swage2019(4,2:6),2)+sum(Swage2019(6,2:6),2)) % business receipts less than 1. Mio
                    100.*sum(Swage2019(4,7:8),2)./(sum(Swage2019(4,7:8),2)+sum(Swage2019(6,7:8),2)) % business receipt between 1 and 5 Mio
                    100.*sum(Swage2019(4,9),2)./(sum(Swage2019(4,9),2)+sum(Swage2019(6,9),2)) % business receipt between 5 and 10 Mio
                    100.*sum(Swage2019(4,10:12),2)./(sum(Swage2019(4,10:12),2)+sum(Swage2019(6,10:12),2))]; % business receipt larger than 10 Mio



histo_Ccorp_2019 = [100.*sum(Cwage2019(4,2:6),2)./(sum(Cwage2019(4,2:6),2)+sum(Cwage2019(6,2:6),2)) % business receipts less than 1. Mio
                    100.*sum(Cwage2019(4,7:8),2)./(sum(Cwage2019(4,7:8),2)+sum(Cwage2019(6,7:8),2)) % business receipt between 1 and 5 Mio
                    100.*sum(Cwage2019(4,9),2)./(sum(Cwage2019(4,9),2)+sum(Cwage2019(6,9),2)) % business receipt between 5 and 10 Mio
                    100.*sum(Cwage2019(4,10:14),2)./(sum(Cwage2019(4,10:14),2)+sum(Cwage2019(6,10:14),2))]; % business receipt larger than 10 Mio


histo_Scorp = [histo_Scorp_2007, histo_Scorp_2013, histo_Scorp_2019];
histo_Ccorp = [histo_Ccorp_2007, histo_Ccorp_2013, histo_Ccorp_2019];


categories = ["< $1 mil." "$1 - 5 mil." "$5 - 10 mil." "> $10 mil."];

fig4 = bar(categories,histo_Scorp);
legend('2007', '2013', '2019')
fig4(1).FaceColor = [0.3 0.3 0.3];
fig4(2).FaceColor = [0.7 0.7 0.7];
fig4(3).FaceColor = [0.9 0.9 0.9];
ylim([0 75])
set(fig4,  'LineWidth',2);
set(gca, 'FontSize',18);
ylabel(' %');
%saveas(gca, 'wageScorpthreeyears','epsc')

%pause;
close all;

fig4 = bar(categories,histo_Ccorp);
%legend('2007', '2013', '2019')
fig4(1).FaceColor = [0.3 0.3 0.3];
fig4(2).FaceColor = [0.7 0.7 0.7];
fig4(3).FaceColor = [0.9 0.9 0.9];
ylim([0 75])
ylabel(' %');
set(fig4,  'LineWidth',2);
set(gca, 'FontSize',18);
%saveas(gca, 'wageCcorpthreeyears','epsc')


%pause;
close all;


fig5 = bar(categories, [histo_Scorp_2013 histo_Ccorp_2013]);
legend('S-corporations', 'C-corporations')
fig5(1).FaceColor = [0.2 0.2 0.2];
fig5(2).FaceColor = [0.6 0.6 0.6];
ylim([0 100])
ylabel(' %');
set(fig5,  'LineWidth',2);
set(gca, 'FontSize',20);
%saveas(gca, 'wage2013','epsc')


%pause;
close all;


% share of businesses in different categories within LFO


business_Scorp_2013 = [100.*sum(Swage2013(7,2:6),2)./Swage2013(7,1) % business receipts less than 1. Mio
                    100.*sum(Swage2013(7,7:8),2)./Swage2013(7,1) % business receipt between 1 and 5 Mio
                    100.*sum(Swage2013(7,9),2)./Swage2013(7,1) % business receipt between 5 and 10 Mio
                    100.*sum(Swage2013(7,10:11),2)./Swage2013(7,1)]; % business receipt larger than 10 Mio



business_Ccorp_2013 = [100.*sum(Cwage2013(7,2:6),2)./Cwage2013(7,1) % business receipts less than 1. Mio
                    100.*sum(Cwage2013(7,7:8),2)./Cwage2013(7,1) % business receipt between 1 and 5 Mio
                    100.*sum(Cwage2013(7,9),2)./Cwage2013(7,1) % business receipt between 5 and 10 Mio
                    100.*sum(Cwage2013(7,10:13),2)./Cwage2013(7,1)]; % business receipt larger than 10 Mio



% share of business receipts in different categories within LFO

receipt_Scorp_2013 = [100.*sum(Swage2013(1,2:6),2)./Swage2013(1,1) % business receipts less than 1. Mio
                    100.*sum(Swage2013(1,7:8),2)./Swage2013(1,1) % business receipt between 1 and 5 Mio
                    100.*sum(Swage2013(1,9),2)./Swage2013(1,1) % business receipt between 5 and 10 Mio
                    100.*sum(Swage2013(1,10:11),2)./Swage2013(1,1)]; % business receipt larger than 10 Mio

receipt_Ccorp_2013 = [100.*sum(Cwage2013(1,2:6),2)./Cwage2013(1,1) % business receipts less than 1. Mio
                    100.*sum(Cwage2013(1,7:8),2)./Cwage2013(1,1) % business receipt between 1 and 5 Mio
                    100.*sum(Cwage2013(1,9),2)./Cwage2013(1,1) % business receipt between 5 and 10 Mio
                    100.*sum(Cwage2013(1,10:13),2)./Cwage2013(1,1)]; % business receipt larger than 10 Mio



fig6 = bar(categories, [business_Scorp_2013 business_Ccorp_2013]);
%legend('S-corporations', 'C-corporations')
fig6(1).FaceColor = [0.2 0.2 0.2];
fig6(2).FaceColor = [0.6 0.6 0.6];
ylim([0 100])
set(fig6,  'LineWidth',2);
set(gca, 'FontSize',20);
ylabel(' %');
%saveas(gca, 'sharebusinessLFO','epsc')

%pause;
close all;

fig7 = bar(categories, [receipt_Scorp_2013 receipt_Ccorp_2013]);
%legend('S-corporation', 'C-corporation')
fig7(1).FaceColor = [0.2 0.2 0.2];
fig7(2).FaceColor = [0.6 0.6 0.6];
ylim([0 100])
set(fig7,  'LineWidth',2);
set(gca, 'FontSize',20);
ylabel(' %');
%saveas(gca, 'sharereceiptLFO','epsc')


%pause;
close all;

