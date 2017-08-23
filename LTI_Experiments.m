%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     author: Mario Senden (mario.senden@maastrichtuniversity.nl)     %%%

% This is a reference implementation of the simulation experiments
% described in (Lange, Senden, Radermacher, De Weerd, submitted)
clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             settings                                %%%

OD_0        =   7.5;        % initial orientation difference
Sessions    =   8;          % number of sessions
Reps        =  25;          % number of times each experiment is repeated
Trials      = 480;          % number of trials per session
P           =   0.;         % separability (proportion)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             parameters                              %%%

N           = 512;          % number of neurons
alpha       =  10;          % width of feedforward bias
sigma_ff    =  45;          % gain of spike encoder
J_ff        =    .5;        % forward connection strength
J_rec       =    1;         % recurrent connection strength
a_e         =    2.2;       % exponent exc. connections
a_i         =    1.4;       % exponent inh. connections
c_e         =    1.2025e-3; % normalization exc. connection
c_i         =    1.6875e-3; % normalization inh. connection
k           =    1.47;      % scaling of variance
C           =        .53;   % decision criterion
eta         =    1.5e-9;    % learning rate
t_sim       =     .5;       % simulation time (seconds)
tau         =    1.5e-2;    % membrane time constant (seconds)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                 setup                               %%%

for i=1:3
   Q{i}         = RM(...    % create a model for each of three quadrants 
                    N,...   % (i.e. experiments)
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
                    t_sim,...
                    tau,...
                    OD_0);  
   Exp{i}.Ab    = zeros(Reps,Sessions);
   Exp{i}.At    = zeros(Reps,Sessions); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             experiments                             %%%

% Exp1  (blue:  135°    ->      105° & 165°      -> 135°)
% Exp2  (red:    //     ->      105° & 165°      -> 135°)
% Exp3  (green: 135°    ->          45°          -> 135°)
% 
% Separability is implemented as follows:
% In experiment 1, a proportion (P) of connection changes as a consequence
% of training in part 1 (at 135°) will be fixed and thus not undergo
% changes in response to training at (potentially) interfering
% orientations in part 2.
% 
% In experiment 2, a proportion (P) of connection weights will be fixed
% prior to part 2 such that they are not recruited by training at 
% (potentially) interfering orientations. These weights will subsequently 
% become malleable prior to part 3. These weights were thus 'reserved' for
% training at 135°.


for r=1:Reps
    
    % part 1 (135° - baseline)
    Q{1}.Phi    = 135;
    Q{3}.Phi    = 135;
    for s=1:Sessions
        Q{1}.session();
        Q{3}.session();
        Exp{1}.Ab(r,s)    = Q{1}.mean_JND;
        Exp{3}.Ab(r,s)    = Q{3}.mean_JND;
    end
    
    % part 2a (105° & 45° - interference)
    Q{1}.Phi    = 105;
    Q{2}.Phi    = 105;
    Q{3}.Phi    =  45;
    Q{1}.set_OD();
    Q{2}.set_OD();
    Q{3}.set_OD();
    Q{1}.fix(P);     
    Q{2}.fix(P);
    for s=1:Sessions
        Q{1}.session();
        Q{2}.session();
        Q{3}.session();
    end
    
    % part 2b (165° - interference)
    Q{1}.Phi    = 165;
    Q{2}.Phi    = 165;
    Q{1}.set_OD();
    Q{2}.set_OD();
    for s=1:Sessions
        Q{1}.session();
        Q{2}.session();
    end
    
    % part 3 (135° - test)
    Q{1}.Phi    = 135;
    Q{2}.Phi    = 135;
    Q{3}.Phi    = 135;
    Q{1}.set_OD(Exp{1}.Ab(r,end));
    Q{2}.set_OD();
    Q{3}.set_OD(Exp{3}.Ab(r,end));
    Q{2}.fix(0);
    for s=1:Sessions
        Q{1}.session();
        Q{2}.session();
        Q{3}.session();
        Exp{1}.At(r,s)    = Q{1}.mean_JND;
        Exp{2}.At(r,s)    = Q{2}.mean_JND;
        Exp{3}.At(r,s)    = Q{3}.mean_JND;
    end
    
    Q{1}.reset_weights();
    Q{2}.reset_weights();
    Q{3}.reset_weights();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             plotting                                %%%

Name        = sprintf('separability = %d%%',P*100);
Pos         = [200 200  950 350];

figure('Color','w','Position' ,Pos,'name',Name)

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
