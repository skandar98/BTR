%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     author: Mario Senden (mario.senden@maastrichtuniversity.nl)     %%%

% This is an implementation of the simulation experiments described in
% Lange, G., Senden, M., Radermacher, A., & De Weerd, P. (2020). 
% Interfering with a memory without erasing its trace. 
% Neural Networks, 121, 339–355. 
% https://doi.org/10.1016/j.neunet.2019.09.027

clear all;close all;clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             settings                                %%%

OD_0     =   4.5;           % initial orientation difference
Sessions =   8;             % number of sessions
Reps     =  1;             % number of times each experiment is repeated
Trials   = 100;             % number of trials per session


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
k        =   1.05;          % scaling of variance
C        =   0.53;          % decision criterion
eta      =   1.4e-11;       % learning rate
mu       =   0;             % exponent of power law weight dependence
t_sim    =   0.5;           % simulation time (seconds)
tau      =   1.5e-2;        % membrane time constant (seconds)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                 setup                               %%%

Q = cell(3,1);
Exp = cell(3,1);    Int = cell(3, 1);
qp = cell(3,1);     qr = cell(3,1);     rp = cell(3,1);     rr = cell(3,1);
OD = cell(3, 1);    Phi = cell(3, 1);
for i=1:3
    Q{i} = RM_GN(...   % create a model for each of three quadrants
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
    Exp{i}.Ab = zeros(Reps,Sessions);
    Exp{i}.At = zeros(Reps,Sessions);
    Int{i}.left = zeros(Reps,Sessions);
    Int{i}.right = zeros(Reps,Sessions);
    % q probe
    qp{i}.Ab        = NaN(Reps, Sessions, Trials);
    qp{i}.At        = NaN(Reps, Sessions, Trials);
    qp{i}.left      = NaN(Reps, Sessions, Trials);
    qp{i}.right     = NaN(Reps, Sessions, Trials);
    % q reference
    qr{i}.Ab        = NaN(Reps, Sessions, Trials);
    qr{i}.At        = NaN(Reps, Sessions, Trials);
    qr{i}.left      = NaN(Reps, Sessions, Trials);
    qr{i}.right     = NaN(Reps, Sessions, Trials);
    % activation probe at end of trial
    rp{i}.Ab        = NaN(Reps, Sessions, Trials, N);
    rp{i}.At        = NaN(Reps, Sessions, Trials, N);
    rp{i}.left      = NaN(Reps, Sessions, Trials, N);
    rp{i}.right     = NaN(Reps, Sessions, Trials, N);
    % activation reference at end of trial
    rr{i}.Ab        = NaN(Reps, Sessions, Trials, N);
    rr{i}.At        = NaN(Reps, Sessions, Trials, N);
    rr{i}.left      = NaN(Reps, Sessions, Trials, N);
    rr{i}.right     = NaN(Reps, Sessions, Trials, N);
    % orientation difference at each trial
    OD{i}.Ab        = NaN(Reps, Sessions, Trials);
    OD{i}.At        = NaN(Reps, Sessions, Trials);
    OD{i}.left      = NaN(Reps, Sessions, Trials);
    OD{i}.right     = NaN(Reps, Sessions, Trials);

    % orientation difference at each trial
    Phi{i}.Ab        = NaN(Reps, Sessions, Trials);
    Phi{i}.At        = NaN(Reps, Sessions, Trials);
    Phi{i}.left      = NaN(Reps, Sessions, Trials);
    Phi{i}.right     = NaN(Reps, Sessions, Trials);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             experiments                             %%%

% Exp1  (green: 135°    ->          45°          -> 135°)
% Exp2  (blue:  135°    ->      105° & 165°      -> 135°)
% Exp3  (red:    //     ->      105° & 165°      -> 135°)
rate_curve = zeros(Reps,5,181,N);
for r=1:Reps
    fprintf('\n - participant %.2d',r)
    % part 1 (135° - baseline)
%     Q{1}.set_PHI(135);
    for i=-90:90
       [rate_curve(r,1,i+91,:), v] = Q{2}.get_response(i); %  tuning curve
    end  
    Q{2}.set_PHI(135);
    for s=1:Sessions
        [qp{2}.Ab(r, s, :), qr{2}.Ab(r, s, :), rp{2}.Ab(r, s, :, :), rr{2}.Ab(r, s, :, :), OD{2}.Ab(r, s, :), Phi{2}.Ab(r, s, :)] = Q{2}.session();
%         Exp{1}.Ab(r,s) = Q{1}.get_JND * 2;
        Exp{2}.Ab(r,s) = Q{2}.get_JND * 2;
    end
    for i=-90:90
       [rate_curve(r,2,i+91,:), v] = Q{2}.get_response(i); %  tuning curve
    end  
    % part 2a (105° & 45° - interference)
%     Q{1}.set_PHI(45);
    Q{2}.set_PHI(105);
%     Q{3}.set_PHI(105);
%     Q{1}.set_OD();
    Q{2}.set_OD();
%     Q{3}.set_OD();
    for s=1:Sessions
%         Q{1}.session();
        [qp{2}.left(r, s, :), qr{2}.left(r, s, :), rp{2}.left(r, s, :, :), rr{2}.left(r, s, :, :), OD{2}.left(r, s, :), Phi{2}.left(r, s, :)] = Q{2}.session();
%         Q{3}.session();
%         Int{1}.left(r,s) = Q{1}.get_JND * 2;
        Int{2}.left(r,s) = Q{2}.get_JND * 2;
%         Int{3}.left(r,s) = Q{3}.get_JND * 2;
    end
    for i=-90:90
       [rate_curve(r,3,i+91,:), v] = Q{2}.get_response(i); %  tuning curve
    end  
%     % part 2b (165° - interference)
    Q{2}.set_PHI(165);
%     Q{3}.set_PHI(165);
    Q{2}.set_OD();
%     Q{3}.set_OD();
    for s=1:Sessions
        [qp{2}.right(r, s, :), qr{2}.right(r, s, :), rp{2}.right(r, s, :, :), rr{2}.right(r, s, :, :), OD{2}.right(r, s, :), Phi{2}.right(r, s, :)] = Q{2}.session();
%         Q{3}.session();
        Int{2}.right(r,s) = Q{2}.get_JND * 2;
%         Int{3}.right(r,s) = Q{3}.get_JND * 2;
%         
    end
    for i=-90:90
       [rate_curve(r,4,i+91,:), v] = Q{2}.get_response(i); %  tuning curve
    end  
    % part 3 (135° - test)
%     Q{1}.set_PHI(135);
    Q{2}.set_PHI(135);
%     Q{3}.set_PHI(135);
%     Q{1}.set_OD(Exp{1}.Ab(r,end));
    Q{2}.set_OD(Exp{2}.Ab(r,end));
%     Q{3}.set_OD();
    for s=1:Sessions
%         Q{1}.session();
        [qp{2}.At(r, s, :), qr{2}.At(r, s, :), rp{2}.At(r, s, :, :), rr{2}.At(r, s, :, :), OD{2}.At(r, s, :), Phi{2}.At(r, s, :)] = Q{2}.session();
%         Q{3}.session();
%         Exp{1}.At(r,s) = Q{1}.get_JND * 2;
        Exp{2}.At(r,s) = Q{2}.get_JND * 2;
%         Exp{3}.At(r,s) = Q{3}.get_JND * 2;
    end
    for i=-90:90
       [rate_curve(r,5,i+91,:), v] = Q{2}.get_response(i); %  tuning curve
    end   
%     Q{1}.reset();
    Q{2}.reset();
%     Q{3}.reset();
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                             plotting                                %%%

% Pos = [200 200  950 350];
% figure('Color','w','Position' ,Pos)

% experiment 1
% subplot(1,3,1,'Fontsize',9)
% hold all
% figure(1)
% hold all
% plot(mean(Exp{1}.Ab),'color',[0 .75 0],'linestyle','--','linewidth',2.5)
% plot(mean(Exp{1}.At),'color',[0 .75 0],'linewidth',2.5)
% set(gca, 'XTick', 1:8)
% set(gca, 'YScale', 'log')
% xlim([0.5 8.5])
% ylim([1.5 8.5])
% xlabel('session')
% title('Experiment 1 (ACA)')
% legend('A_B','A_T')
% legend('boxoff')

% experiment 2
figure(2) % JND reference
subplot(121)
hold all, grid on
plot(mean(Exp{2}.Ab),'color',[0 0 .75],'linestyle','--','linewidth',2.5)
plot(mean(Exp{2}.At),'color',[0 0 .75],'linewidth',2.5)
set(gca, 'XTick', 1:8), set(gca, 'YScale', 'log')
xlim([0.5 8.5]),        ylim([1.5 8.5])
xlabel('session'),      ylabel('JND [deg]')
title('Experiment 2 (ABA)'),    legend('A_B','A_T'),        legend('boxoff')

subplot(122) % interference
hold all, grid on 
plot(mean(Int{2}.left),'r-.', 'linewidth',2.5)
plot(mean(Int{2}.right),'g--', 'linewidth',2.5)
set(gca, 'XTick', 1:8),     set(gca, 'YScale', 'log')
xlim([0.5 8.5]),            ylim([1.5 8.5])
xlabel('session'),          ylabel('JND [deg]')
title('Experiment 2 (ABA)')
legend('$105^\circ$','$165^\circ$', 'interpreter', 'latex'),    legend('boxoff')

figure(3) % q(sessions) probe
subplot(121)
plot(mean(squeeze(qp{2}.Ab(1, :, :)), 2), 'color',[0 0 .75],'linestyle','--','linewidth',2.5), hold on
plot(mean(squeeze(qp{2}.At(1, :, :)), 2), 'color',[0 0 .75],'linewidth',2.5), grid on, hold off
xlabel('session'), ylabel('q [-]')
xlim([0.5 8.5])

subplot(122)
plot(mean(squeeze(qp{2}.left(1, :, :)), 2), 'r-.', 'linewidth',2.5), grid on, hold on
plot(mean(squeeze(qp{2}.right(1, :, :)), 2),'g--', 'linewidth',2.5), hold off
xlabel('session'), ylabel('q [-]')
xlim([0.5 8.5])
legend('$105^\circ$','$165^\circ$', 'interpreter', 'latex')

figure(4) % q(sessions) reference
subplot(121)
plot(mean(squeeze(qr{2}.Ab(1, :, :)), 2), 'color',[0 0 .75],'linestyle','--','linewidth',2.5), hold on 
plot(mean(squeeze(qr{2}.At(1, :, :)), 2), 'color',[0 0 .75],'linewidth',2.5), grid on
xlabel('session'), ylabel('q [-]')
xlim([0.5 8.5])

subplot(122)
plot(mean(squeeze(qr{2}.left(1, :, :)), 2), 'r-.', 'linewidth',2.5), grid on, hold on
plot(mean(squeeze(qr{2}.right(1, :, :)), 2),'g--', 'linewidth',2.5), hold off
xlabel('session'), ylabel('q [-]')
xlim([0.5 8.5])
legend('$105^\circ$','$165^\circ$', 'interpreter', 'latex')

%% Part A
symbols = ['+', '*', 's', '^', 'd', 'p', 'h', '<'];     % symbols per depicted session (max 8)
legend_entries = ["Session1", "Session2", "Session3", "Session4", "Session5", "Session6", "Session7", "Session8"]';
legend_indices = [1 3 5 7 9 11 13 15];
% PCA of reference dataset
RR = reshape(squeeze(rr{2}.Ab(1, :, :, :)), Sessions*Trials, N);
[coeff_rr,score_rr,latent_rr, ~, explained_rr, mu_rr] = pca(RR);
% explained variance
figure(51)
    plot(cumsum(explained_rr(1:5)), 'LineWidth', 1.5), grid on
    xlabel('# PC'), ylabel('CEV [%]')

figure(5) % colored by kinetic energy
k = 1;  hp = NaN(Sessions, 1); hold all
for i = 1:Sessions
    % Reference
    score = (squeeze(rr{2}.Ab(1, i, :, :)) - mu_rr)*coeff_rr;
    hp(2*i) = scatter(score(:, 2), score(:, 1), [], squeeze(qr{2}.Ab(1, i, :)), 'Marker', '.');

    % Probe
    score = (squeeze(rp{2}.Ab(1, i, :, :)) - mu_rr)*coeff_rr;
    hp(2*i-1) = scatter(score(:, 2), score(:, 1), [], squeeze(qp{2}.Ab(1, i, :)), 'Marker',symbols(k));
    k = k+1;
end
xlabel('PC 2'), ylabel('PC 1'), grid on, colorbar(), hold off
legend(hp(legend_indices), legend_entries)

figure(52) % colored by orientation difference
k = 1;  hp = NaN(Sessions, 1); hold all
for i = 1:Sessions
    % Reference
    score = (squeeze(rr{2}.Ab(1, i, :, :)) - mu_rr)*coeff_rr;
    hp(2*i) = scatter(score(:, 2), score(:, 1), 'k.');

    % Probe
    score = (squeeze(rp{2}.Ab(1, i, :, :)) - mu_rr)*coeff_rr;
    hp(2*i-1)= scatter(score(:, 2), score(:, 1), [], squeeze(Phi{2}.Ab(1, i, :)), 'Marker',symbols(k)); % potentially also OD
    k = k+1; 
end
xlabel('PC 2'), ylabel('PC 1'), grid on, colorbar(), hold off
legend(hp(legend_indices), legend_entries)



%% Part B
% interference left 105
RRl = reshape(squeeze(rr{2}.left(1, :, :, :)), Sessions*Trials, N);
[coeff_rrl,score_rrl,latent_rrl, ~, explained_rrl, mu_rrl] = pca(RRl);
% explained variance
figure(71)
plot(cumsum(explained_rrl(1:5)), 'LineWidth', 1.5), grid on
xlabel('# PC'), ylabel('CEV [%]')

figure(7) % colored by kinetic energy
k = 1;  hp = NaN(Sessions, 1); hold all
for i = 1:Sessions
    % Reference
    score = (squeeze(rr{2}.left(1, i, :, :)) - mu_rrl)*coeff_rrl;
    hp(2*i) = scatter(score(:, 2), score(:, 1),  [], squeeze(qr{2}.left(1, i, :)), 'Marker', '.');

    % Probe
    score = (squeeze(rp{2}.left(1, i, :, :)) - mu_rrl)*coeff_rrl;
    hp(2*i-1)= scatter(score(:, 2), score(:, 1), [], squeeze(qp{2}.left(1, i, :)), 'Marker',symbols(k));
    k = k+1;
end
xlabel('PC 2'), ylabel('PC 1'), grid on, colorbar(), hold off
legend(hp(legend_indices), legend_entries)

figure(72) % colored by orientation difference
k = 1;  hp = NaN(Sessions, 1); hold all
for i = 1:Sessions
    % Reference
    score = (squeeze(rr{2}.left(1, i, :, :)) - mu_rrl)*coeff_rrl;
    hp(2*i) = scatter(score(:, 2), score(:, 1), 'k.');

    % Probe
    score = (squeeze(rp{2}.left(1, i, :, :)) - mu_rrl)*coeff_rrl;
    hp(2*i-1)= scatter(score(:, 2), score(:, 1), [], squeeze(Phi{2}.left(1, i, :)), 'Marker',symbols(k)); % potentially also OD
    k = k+1;
end
xlabel('PC 2'), ylabel('PC 1'), grid on, colorbar(), hold off
legend(hp(legend_indices), legend_entries)

% interference right 165
RRr = reshape(squeeze(rr{2}.right(1, :, :, :)), Sessions*Trials, N);
[coeff_rrr,score_rrr,latent_rrr, ~, explained_rrr, mu_rrr] = pca(RRr);
% explained variance
figure(81)
plot(cumsum(explained_rrr(1:5)), 'LineWidth', 1.5), grid on
xlabel('# PC'), ylabel('CEV [%]')

figure(8) % colored by kinetic energy
k = 1;  hp = NaN(Sessions, 1); hold all
for i = 1:Sessions
    % Reference
    score = (squeeze(rr{2}.right(1, i, :, :)) - mu_rrr)*coeff_rrr;
    hp(2*i) = scatter(score(:, 2), score(:, 1), [], squeeze(qr{2}.right(1, i, :)), 'Marker', '.');

    % Probe
    score = (squeeze(rp{2}.right(1, i, :, :)) - mu_rrr)*coeff_rrr;
    hp(2*i-1)= scatter(score(:, 2), score(:, 1), [], squeeze(qp{2}.right(1, i, :)), 'Marker',symbols(k));
    k = k+1;
end
xlabel('PC 2'), ylabel('PC 1'), grid on, colorbar(), hold off
legend(hp(legend_indices), legend_entries)


figure(82) % colored by orientation difference
k = 1;  hp = NaN(Sessions, 1); hold all
for i = 1:Sessions
    % Reference
    score = (squeeze(rr{2}.right(1, i, :, :)) - mu_rrr)*coeff_rrr;
    hp(2*i) = scatter(score(:, 2), score(:, 1), 'k.'); 

    % Probe
    score = (squeeze(rp{2}.right(1, i, :, :)) - mu_rrr)*coeff_rrr;
    hp(2*i-1)= scatter(score(:, 2), score(:, 1), [], squeeze(Phi{2}.right(1, i, :)), 'Marker',symbols(k)); % potentially also OD
    k = k+1;
end
xlabel('PC 2'), ylabel('PC 1'), grid on, colorbar(), hold off
legend(hp(legend_indices), legend_entries)

%% Part A (At)
RRt = reshape(squeeze(rr{2}.At(1, :, :, :)), Sessions*Trials, N);
[coeff_rrt,score_rrt,latent_rrt, ~, explained_rrt, mu_rrt] = pca(RRt);
% explained variance
figure(91)
plot(cumsum(explained_rrt(1:5)), 'LineWidth', 1.5), grid on
xlabel('# PC'), ylabel('CEV [%]')

figure(9) % colored by kinetic energy
k = 1;  hp = NaN(Sessions, 1); hold all
for i = 1:Sessions
    % Reference
    score = (squeeze(rr{2}.At(1, i, :, :)) - mu_rrt)*coeff_rrt;
    hp(2*i) = scatter(score(:, 2), score(:, 1), [], squeeze(qr{2}.At(1, i, :)), 'Marker', '.'); 

    % Probe
    score = (squeeze(rp{2}.At(1, i, :, :)) - mu_rrt)*coeff_rrt;
    hp(2*i-1)= scatter(score(:, 2), score(:, 1), [], squeeze(qp{2}.At(1, i, :)), 'Marker',symbols(k));
    k = k+1;
end
xlabel('PC 2'), ylabel('PC 1'), grid on, colorbar(), hold off
legend(hp(legend_indices), legend_entries)


figure(92) % colored by orientation difference
k = 1;  hp = NaN(Sessions, 1); hold all
for i = 1:Sessions
    % Reference
    score = (squeeze(rr{2}.At(1, i, :, :)) - mu_rrt)*coeff_rrt;
    hp(2*i) = scatter(score(:, 2), score(:, 1), 'k.'); 

    % Probe
    score = (squeeze(rp{2}.At(1, i, :, :)) - mu_rrt)*coeff_rrt;
    hp(2*i-1)= scatter(score(:, 2), score(:, 1), [], squeeze(Phi{2}.At(1, i, :)), 'Marker',symbols(k)); % potentially also OD
    k = k+1;
end
xlabel('PC 2'), ylabel('PC 1'), grid on, colorbar(), hold off
legend(hp(legend_indices), legend_entries)

%% TC plot
figure(1) 
X   = linspace(45, 225, 181);
%X   = linspace(45,225,100) 
idx = round(linspace(40,472,7));


subplot(3,3,2)
for i = 1:7
    plot(X,squeeze(rate_curve(1,1,:,idx(i))))
    hold all
end
xlim([45 225])
title("naive")
xlabel("Orientation (degrees)")
ylabel("Firing rate")

subplot(3,3,3)
for i = 1:7
    plot(X,squeeze(rate_curve(1,2,:,idx(i))))
    hold all
end
xlim([45 225])
xline(135,'r')
title("AB (135)")
xlabel("Orientation (degrees)")
ylabel("Firing rate")

subplot(3,3,4)
for i = 1:7
    plot(X,squeeze(rate_curve(1,3,:,idx(i))))
    hold all
end
xlim([45 225])
xline(105,'r')
title("Left (105)")
xlabel("Orientation (degrees)")
ylabel("Firing rate")

subplot(3,3,6)
for i = 1:7
    plot(X,squeeze(rate_curve(1,4,:,idx(i))))
    hold all
end
xlim([45 225])
xline(165,'r')
title("Right (165)")
xlabel("Orientation (degrees)")
ylabel("Firing rate")

subplot(3,3,8)
for i = 1:7
    plot(X,squeeze(rate_curve(1,5,:,idx(i))))
    hold all
end
xlim([45 225])
xline(135,'r')
title("AT (135)")
xlabel("Orientation (degrees)")
ylabel("Firing rate")




