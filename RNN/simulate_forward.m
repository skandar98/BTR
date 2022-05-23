clc; clear; close all;

%declare input and output path
path_in ="\BTR\FPAnalysis\simulation_results\";
path_out = "\BTR\FPAnalysis\analysis_results\"; 

%load input variables
load(path_in+'inputs1.mat');
load(path_in+'fp.mat')
fp = fp';
rec1 = V_rec(:,4);
ff1 = V_ff(:,4);
W1 = W(:,:,4);

t_sim = .5;

dV = @(V_rec,V_ff,W,alpha,tau)...
      (-V_rec+V_ff+W*alpha*max(V_rec,0))/tau;


[~,v] = ode45(@(t,v)dV(v, ff1, ...
            W1, alpha, tau),...
           [0 t_sim], fp);

f=figure;
mesh(v);
saveas(f,'simulate.png');
