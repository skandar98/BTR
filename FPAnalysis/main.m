clc; clear; close all;

%declare input and output path
path_in ="C:\Users\skand\OneDrive\Dokumentumok\MATLAB\BTR\FPAnalysis\simulation_results\";
path_out = "C:\Users\skand\OneDrive\Dokumentumok\MATLAB\BTR\FPAnalysis\analysis_results\"; 

%load input variables
load(path_in+'inputs1.mat');

stepSize = 1e-3;
beta1 = 0.9;
beta2 = 0.99;
epsilon = 1e-7;

%
W1 = W(:,:,4);
V_ff1 = V_ff(:,4);


x = dlarray(V_rec(:,4));


[y, grad] = dlfeval(@gradFun,x,V_ff1, W1, alpha, tau);

%[x1, fval, exitflag, output] = fmin_adam(dlfeval(@gradFun,x,V_ff1, W1, alpha, tau), x);

function y = qfun(V_rec, V_ff1, W1, alpha, tau)
    dV = (-V_rec+V_ff1+W1.*alpha.*max(V_rec,0))./tau;
    y =sum(dV.^2, 'all')./512;
end

function [y,grad] = gradFun(x, V_ff1, W1, alpha, tau)
    y = qfun(x, V_ff1, W1, alpha, tau);
    grad = dlgradient(y,x);
end

