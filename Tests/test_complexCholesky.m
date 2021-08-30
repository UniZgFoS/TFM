clc;
close all;
clear;
addpath('../Utility/');

n = 30;
A = randn(n,n) + 1j*randn(n,n);
A = 0.5 * ( A + transp(A) );


R = function_complexCholesky( A );
norm( R*transp(R) - A ) % A = R * R^T
norm( triu(R,+1) ) % R is lower trinagular




