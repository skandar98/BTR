clear;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               settings                                 %

path ="\BTR\FPAnalysis\simulation_results\";


OD_0     =   4.5;           % initial orientation difference
Sessions =   4;             % number of sessions
Reps     =   1;             % number of times each experiment is repeated
Trials   = 50;              % number of trials per session

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              parameters                                %

N        = 512;             % number of neurons
alpha    =  10;             % gain of spike encoder
sigma_ff =  45;             % width of feedforward bias
J_ff     =   0.5;           % forward connection strength
J_rec    =   1;             % recurrent connection strength
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



Q= RM(...   % create a model for each of three quadrants
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
W = zeros(512, 512, Sessions*3);
rate = zeros(512, Sessions*3);
V_rec = zeros(512, Sessions*3);
V_ff = zeros(512, Sessions*3);
tic
for r=1:Reps
    Q.set_PHI(135);
    for s=1:Sessions
        W(:,:,s) = Q.session();
        [V_ff(:,s), V_rec(:,s), rate(:,s)]= Q.get_response();
        Exp.Ab(r, s) = Q.get_JND * 2;

    end
    
     Q.set_PHI(45);
     Q.set_OD();
     
     for s=1:Sessions
          W(:,:,s+Sessions) = Q.session();
          [V_ff(:,s+Sessions), V_rec(:,s+Sessions), rate(:,s)]= Q.get_response();
          Int(r, s) = Q.get_JND * 2;
          V_ff(:,s+Sessions) 
     end
%     
%     Q.set_PHI(135);
%     Q.set_OD();
%     for s=1:Sessions
%         W(s) = Q.session();
%         Exp.At(r, s) = Q.get_JND * 2;
%     end
%     
    Q.reset();
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                save                                   %%

save(path+'inputs1.mat', 'W', 'rate', 'V_rec', "V_ff",'alpha','tau');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                              plotting                                 %%

% Pos = [200 200  950 350];
figure('Color','w')

% experiment 1
% subplot(1,3,1,'Fontsize',9)
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
