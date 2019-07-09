%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     author: Mario Senden (mario.senden@maastrichtuniversity.nl)     %%%

% This is an implementation of the simulation experiments
% described in: Lange G, Senden M, Radermacher A, De Weerd P.
% Interfering with a memory without erasing its trace (submitted).
clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             settings                                %%%

OD_0     =   7.5;           % initial orientation difference
Sessions =   8;             % number of sessions
Reps     =   4;             % number of times each experiment is repeated
Trials   = 480;             % number of trials per session
P        = 0.75;


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

for i=1:4
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
    Q{i}.set_PHI(135);
    Exp{i}      = zeros(Reps,Sessions);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             experiments                             %%%

% Exp1  (correct = 0; incorrect = 0)
% Exp2  (correct = P; incorrect = 0)
% Exp3  (correct = 0; incorrect = P)
% Exp4  (correct = P; incorrect = P)

for r=1:Reps
    fprintf('\n - participant %.2d',r)
    
    Q{2}.set_CORRECT(P);
    Q{3}.set_INCORRECT(P);
    Q{4}.set_CORRECT(P);
    Q{4}.set_INCORRECT(P);
    for s=1:Sessions
        for i=1:4
            Q{i}.session();
            Exp{i}(r,s)  = Q{i}.get_JND;
        end
    end
    
    
    
    for i=1:4
        [We_trained{i},Wi_trained{i}] = Q{i}.get_weights();
        Q{i}.reset();
        [We{i},Wi{i}] = Q{i}.get_weights();
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             plotting                                %%%

Pos = [200 200  950 350];
figure('Color','w','Position' ,Pos)

% experiment 1
hold all
plot(mean(Exp{1}),'color',[0 0 .75],'linewidth',2.5)
plot(mean(Exp{2}),'color',[.75 0 0],'linewidth',2.5)
plot(mean(Exp{3}),'color',[0 .75 0],'linewidth',2.5)
plot(mean(Exp{4}),'color',[.75 .75 0],'linewidth',2.5)
set(gca, 'XTick', 1:8)
set(gca, 'YScale', 'log')
xlim([0.5 8.5])
ylim([1.5 8.5])
xlabel('session')
ylabel('JND (degree)')
title('probabilistic feedback on')
legend('neither trials','incorrect trials', 'correct trials', 'both trials')
legend('boxoff')

%%

