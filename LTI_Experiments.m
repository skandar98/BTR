%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     author: Mario Senden (mario.senden@maastrichtuniversity.nl)     %%%

clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             settings                                %%%

Experiment  =   1;
reserved    =   0;
JND_0       =   7.5;
Subjects    =  25;
Sessions    =   8;
Staircases  =   4;
Trials      = 120;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             parameters                              %%%

N           = 512;                  % number of neurons
sigma_ff    =  45;                  % width of feedforward orientation bias
alpha       =  10;                  % gain of spike encoder
J_ff        =    .5;                % feedforward connection strength
J_rec       =   1;                  % recurrent connection strength
a_e         =   2.2;                % exponent excitatory connections
a_i         =   1.4;                % exponent inhibitory connections
c_e         =    .00118;            % normalization excitatory connections
c_i         =    .00165;            % normalization inhibitory connections
k           =   1.47;               % scaling of variance
C           =    .53;               % decision criterion
eta         = 15e-9;                % learning rate
t_sim       =    .5;                % simulation time
dt          =  2e-3;                % time step in numerical integration
tau         = 15e-3;                % membrane time constant


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                               setup                                 %%%

JNDs        = ones(Subjects,...     %
            Sessions)*JND_0;
V1          = RM(N,...
                 sigma_ff,...
                 alpha,...
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
                 dt,...
                 tau);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             experiments                             %%%
fprintf('> Experiment: %d\n', Experiment)
switch Experiment
    case 1
        fprintf('> Color code: blue \n')
        fprintf('> Sequence: 0° -> %s30° -> 0° \n',177)
        fprintf('> Connections reserved for learning around 0°: %d%s \n',...
            reserved*100,37)
        for subj=1:Subjects
            
            for sess=1:Sessions
                for t=Trials*Staircases
                    
                end
            end
        end
        
    case 2
        fprintf('> Color code: blue \n')
        fprintf('> Sequence: %s30° -> 0° \n',177)
        fprintf('> Connections reserved for learning around 0°: %d%s \n',...
            reserved*100,37)
        for subj=1:Subjects
            
            for sess=1:Sessions
                for t=Trials*Staircases
                    
                end
            end
        end
        
    case 3
        fprintf('> Color code: blue \n')
        fprintf('> Sequence: 0° -> -90° -> 0° \n')
        fprintf('> Connections reserved for learning around 0°: %d%s \n',...
            reserved*100,37)
        for subj=1:Subjects
            
            for sess=1:Sessions
                for t=Trials*Staircases
                    
                end
            end
        end

end
V1.Phi=30;
V1.trial()
