clear; clc; close all
format long g

n = 100;
a_min = 1e-6;
a_max = 50;

% Log grid

a_grid = linspace(log(a_min),log(a_max),n)';
a_grid = exp(a_grid);

figure
plot(a_grid,a_grid,'-o')
title('Log-spaced grid')

% Use my routine
a_space = 10;
a_grid1 = make_grid(a_min,a_max,n,a_space,1);

figure
plot(a_grid1,a_grid1,'-o')
title('Left-skewed grid')






