%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     author: Mario Senden (mario.senden@maastrichtuniversity.nl)     %%%

% This is an implementation of the simulation experiments
% described in: Lange G, Senden M, Radermacher A, De Weerd P.
% Interfering with a memory without disrupting its trace (submitted).
clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             settings                                %%%

OD_0     =   7.5;           % initial orientation difference
Sessions =   8;             % number of sessions
Reps     =  25;             % number of times each experiment is repeated
Trials   = 480;             % number of trials per session


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             parameters                              %%%

N        = 512;             % number of neurons
alpha    =  10;             % gain of spike encoder
sigma_ff =  45;             % width of feedforward bias
J_ff     =   0.5;           % forward connection strength
J_rec    =   1;             % recurrent connection strength
a_e      =   2.2;           % exponent exc. connections
a_i      =   1.4;           % exponent inh. connections
c_e      =   1.2025e-3;     % normalization exc. connection
c_i      =   1.6875e-3;     % normalization inh. connection
k        =   4;             % scaling of variance
C        =   0.53;          % decision criterion
eta      =   1.4e-11;       % learning rate
mu       =   0;             % exponent of power law weight dependence
t_sim    =   0.5;           % simulation time (seconds)
tau      =   1.5e-2;        % membrane time constant (seconds)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                 setup                               %%%

for i=1:3
    Q{i}         = RM(...   % create a model for each of three quadrants
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
    Exp{i}.Ab    = zeros(Reps,Sessions);
    Exp{i}.At    = zeros(Reps,Sessions);
    Int{i}.left  = zeros(Reps,Sessions);
    Int{i}.right = zeros(Reps,Sessions);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             experiments                             %%%

% Exp1  (blue:  135°    ->      105° & 165°      -> 135°)
% Exp2  (red:    //     ->      105° & 165°      -> 135°)
% Exp3  (green: 135°    ->          45°          -> 135°)

for r=1:Reps
    fprintf('\n - participant %.2d',r)
    % part 1 (135° - baseline)
    Q{1}.set_PHI(135);
    Q{3}.set_PHI(135);
    for s=1:Sessions
        Q{1}.session();
        Q{3}.session();
        Exp{1}.Ab(r,s)  = Q{1}.get_JND;
        Exp{3}.Ab(r,s)  = Q{3}.get_JND;
    end
    
    % part 2a (105° & 45° - interference)
    Q{1}.set_PHI(105);
    Q{2}.set_PHI(105);
    Q{3}.set_PHI(45);
    Q{1}.set_OD();
    Q{2}.set_OD();
    Q{3}.set_OD();
    for s=1:Sessions
        Q{1}.session();
        Q{2}.session();
        Q{3}.session();
        Int{1}.left     = Q{1}.get_JND;
        Int{2}.left     = Q{2}.get_JND;
        Int{3}.left     = Q{3}.get_JND;
    end
    
    % part 2b (165° - interference)
    Q{1}.set_PHI(165);
    Q{2}.set_PHI(165);
    Q{1}.set_OD();
    Q{2}.set_OD();
    for s=1:Sessions
        Q{1}.session();
        Q{2}.session();
        Int{1}.right    = Q{1}.get_JND;
        Int{2}.right    = Q{2}.get_JND;
        
    end
    
    % part 3 (135° - test)
    Q{1}.set_PHI(135);
    Q{2}.set_PHI(135);
    Q{3}.set_PHI(135);
    Q{1}.set_OD(Exp{1}.Ab(r,end));
    Q{2}.set_OD();
    Q{3}.set_OD(Exp{3}.Ab(r,end));
    for s=1:Sessions
        Q{1}.session();
        Q{2}.session();
        Q{3}.session();
        Exp{1}.At(r,s)  = Q{1}.get_JND;
        Exp{2}.At(r,s)  = Q{2}.get_JND;
        Exp{3}.At(r,s)  = Q{3}.get_JND;
    end
    
    Q{1}.reset();
    Q{2}.reset();
    Q{3}.reset();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             plotting                                %%%

Pos = [200 200  950 350];
figure('Color','w','Position' ,Pos)

% experiment 1
subplot(1,3,1,'Fontsize',9)
hold all
plot(mean(Exp{1}.Ab),'color',[0 0 .75],'linestyle','--','linewidth',2.5)
plot(mean(Exp{1}.At),'color',[0 0 .75],'linewidth',2.5)
set(gca, 'XTick', 1:8)
set(gca, 'YScale', 'log')
xlim([0.5 8.5])
ylim([1.5 8.5])
xlabel('session')
ylabel('JND (degree)')
title('Experiment 1 (ABA)')
legend('A_B','A_T')
legend('boxoff')

% experiment 2
subplot(1,3,2,'Fontsize',9)
hold all
plot(mean(Exp{1}.Ab),'color',[0 0 .75],'linestyle','--','linewidth',2.5)
plot(mean(Exp{2}.At),'color',[.75 0 0],'linewidth',2.5)
set(gca, 'XTick', 1:8)
set(gca, 'YScale', 'log')
xlim([0.5 8.5])
ylim([1.5 8.5])
xlabel('session')
title('Experiment 2 (BA)')
legend('A_B','A_T')
legend('boxoff')

% experiment 3
subplot(1,3,3,'Fontsize',9)
hold all
plot(mean(Exp{3}.Ab),'color',[0 .75 0],'linestyle','--','linewidth',2.5)
plot(mean(Exp{3}.At),'color',[0 .75 0],'linewidth',2.5)
set(gca, 'XTick', 1:8)
set(gca, 'YScale', 'log')
xlim([0.5 8.5])
ylim([1.5 8.5])
xlabel('session')
title('Experiment 3 (ACA)')
legend('A_B','A_T')
legend('boxoff')