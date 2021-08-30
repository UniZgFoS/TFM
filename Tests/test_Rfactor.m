clear;
clc;
close all;
addpath('../Utility/');

m = 100;
n = 60;
A = randn(m,n) + 1j*randn(m,n);

R = function_Rfactor(A);

v = randn(n,1) + 1j*randn(n,1);

s1 = sum( (A*v).^2 )
s2 = sum( (R*v).^2 )

abs(s1-s2)/abs(s1)







