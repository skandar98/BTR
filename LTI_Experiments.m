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
Reps     =   4;             % number of times each experiment is repeated
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

% Exp1  (blue:  135°    ->      105° & 165°      -> 135°)
% Exp2  (red:    //     ->      105° & 165°      -> 135°)
% Exp3  (green: 135°    ->          45°          -> 135°)

for r=1:Reps
    fprintf('\n - participant %.2d',r)
    
    Q{2}.set_CORRECT(0.5);
    Q{3}.set_INCORRECT(0.5);
    Q{4}.set_CORRECT(0.5);
    Q{4}.set_INCORRECT(0.5);
    for s=1:Sessions
        for i=1:4
            Q{i}.session();
            Exp{i}(r,s)  = Q{i}.get_JND;
        end
    end
    
    for i=1:4
        Q{i}.reset();
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             plotting                                %%%

Pos = [200 200  950 350];
figure('Color','w','Position' ,Pos)

% experiment 1
hold all
for i=1:4
    plot(mean(Exp{i}),'color',[0 0 .75],'linewidth',2.5)
end
set(gca, 'XTick', 1:8)
set(gca, 'YScale', 'log')
xlim([0.5 8.5])
ylim([1.5 8.5])
xlabel('session')
ylabel('JND (degree)')
title('probabilistic feedback on')
legend('neither trials','incorrect trials', 'correct trials', 'both trials')
legend('boxoff')
