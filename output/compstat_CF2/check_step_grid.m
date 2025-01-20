% compstat_CF2
clear
clc
close all

n = 32;

x_min = 0.396;
x_max = 0.52;
step = (x_max-x_min)/(n-1);

% Now we set x_max = 0.55
% Keep step constant, what is the right n?
x_max = 0.6;
n1 = (x_max-x_min)/step+1;
disp(n1)