clear;
close all;
clc;

addpath('functions/');

x = [-1 0 1];
[z, zx] = legendre_1d(x, 0: 2);
