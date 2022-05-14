clear;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               settings                                 %

path ="C:\Users\skand\OneDrive\Dokumentumok\MATLAB\BTR\FPAnalysis\simulation_results\";


OD_0     =   4.5;           % initial orientation difference
Sessions =   8;             % number of sessions
Reps     =   2;             % number of times each experiment is repeated
Trials   = 450;              % number of trials per session

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              parameters                                %

N        = 512;             % number of neurons
alpha    =  10;             % gain of spike encoder
sigma_ff =  45;             % width of feedforward bias
J_ff     =   0.5;           % forward connection strength
J_rec    =   0.025;          % recurrent connection strength
a_e      =   2.2;           % exponent exc. connections
a_i      =   1.4;           % exponent inh. connections
c_e      =   1.2025e-3;     % normalization exc. connection
c_i      =   1.6875e-3;     % normalization inh. connection
k        =   1.05;          % scaling of variance
C        =   0.53;          % decision criterion
eta      =   1.4e-10;       % learning rate
mu       =   0;             % exponent of power law weight dependence
t_sim    =   0.5;           % simulation time (seconds)
tau      =   1.5e-2;        % membrane time constant (seconds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                setup                                   %



Q= RM_GL(...   % create a model for each of three quadrants
    N,...               % (i.e. experiments)
    alpha,...
    sigma_ff,...
    J_ff,...
    J_rec,...
    a_e,...
    a_i,...
    c_e,...
    c_i,...
    k,...
    C,...
    eta,...
    mu,...
    t_sim,...
    tau,...
    Trials,...
    OD_0);

Exp.Ab = zeros(Reps,Sessions);
Exp.At = zeros(Reps,Sessions);
Int = zeros(Reps, Sessions);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              simulate                                 %%
%W = zeros(512, 512, Trials*Sessions,3);
%V_rec = zeros(512, Trials*Sessions,3);
%V_ff = zeros(512, Trials*Sessions,3);
%C = zeros(Trials*Sessions,3);
rate_curve = NaN(Reps,4,100,N);

tic

for r=1:Reps
    for i=1:100
        [rate_curve(r,1,i,:), v] = Q.get_response(i); %baseline
    end

    Q.set_PHI(135);
    for s=1:Sessions
        idx = (1+Trials*(s-1)):Trials*s;
       % [V_ff(:,idx,1), V_rec(:,idx,1), W(:,:,idx,1),C(idx,1)]= Q.session();
        [V_ff, V_rec, W,C]= Q.session();
        Exp.Ab(r, s) = Q.get_JND * 2;
    end

    for i=1:100
        [rate_curve(r,2,i,:), v] = Q.get_response(i); %Ab tuning curve
    end

     Q.set_PHI(165);
     Q.set_OD();
     
     for s=1:Sessions
          idx = (1+Trials*(s-1)):Trials*s;
          %[V_ff(:,idx,2), V_rec(:,idx,2), W(:,:,idx,2),C(idx,2)]= Q.session();
          [V_ff, V_rec, W,C]= Q.session();
          Int(r, s) = Q.get_JND * 2;
     end
    for i=1:100
        [rate_curve(r,3,i,:), v] = Q.get_response(i); % INT tuning curve
    end    
    Q.set_PHI(135);
    Q.set_OD();
    for s=1:Sessions
        idx = (1+Trials*(s-1)):Trials*s;
        [V_ff, V_rec, W,C]= Q.session();
        %[V_ff(:,idx,3), V_rec(:,idx,3), W(:,:,idx,3),C(idx,3)]= Q.session();
        Exp.At(r, s) = Q.get_JND * 2;
    end
    for i=1:100
        [rate_curve(r,4,i,:), v] = Q.get_response(i); % INT tuning curve
    end    

    Q.reset();
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                save                                   %%
save(path+'inputs_small_W_1rep.mat', 'W', 'V_rec', "V_ff",'alpha','tau','Trials','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              plotting                                 %%

% Pos = [200 200  950 350];
figure('Color','w')

% experiment 1
% subplot(1,3,1,'Fontsize',9)
if Reps == 0
hold all
plot((Exp.Ab),'color',[0 .75 0],'linestyle','--','linewidth',2.5)
plot((Exp.At),'color',[0 .75 0],'linewidth',2.5)
plot(mean(Int), 'color', [0 .75 0], 'linestyle','-.','linewidth', 2.5)
set(gca, 'XTick', 1:Sessions)
set(gca, 'YScale', 'log')
xlim([0.5 Sessions+0.5])
ylim([1.5 8.5])
xlabel('session'), ylabel('JND [deg]')
title('Experiment 1 (ACA)')
legend('A_B','Int','A_T')
legend('boxoff')
else
hold all
plot(mean(Exp.Ab),'color',[0 .75 0],'linestyle','--','linewidth',2.5)
plot(mean(Exp.At),'color',[0 .75 0],'linewidth',2.5)
plot(mean(Int), 'color', [0 .75 0], 'linestyle','-.','linewidth', 2.5)
set(gca, 'XTick', 1:Sessions)
set(gca, 'YScale', 'log')
xlim([0.5 Sessions+0.5])
ylim([1.5 8.5])
xlabel('session'), ylabel('JND [deg]')
title('Experiment 1 (ACA)')
legend('A_B','Int','A_T')
legend('boxoff')
end

figure(2)
X = linspace((135-100), (135+100), 100);
subplot(2,2,1)
for i = 1:7
    plot(X,mean(squeeze(rate_curve(:,1,:,(i-1)*73+1)),1))
    hold all
end
subplot(2,2,2)
for i = 1:7
    plot(X,mean(squeeze(rate_curve(:,2,:,(i-1)*73+1)),1))
    hold all
end
subplot(2,2,3)
for i = 1:7
    plot(X,mean(squeeze(rate_curve(:,3,:,(i-1)*73+1)),1))
    hold all
end
subplot(2,2,4)
for i = 1:7
    plot(X,mean(squeeze(rate_curve(:,4,:,(i-1)*73+1)),1))
    hold all
end


