function [outputArg1,outputArg2] = optimizer(fun,x0)
%Finds fixed points of RNN

path ="C:\Users\skand\OneDrive\Dokumentumok\MATLAB\BTR\FPAnalysis\simulation_results\";

max_iter = 100
lr = 0.001
beta_1 = 0.9
beta_2 = 0.999
epsilon = 1e-3


fun = @(V_rec,V_ff,W,alpha,tau)...
    (-V_rec+V_ff+W*alpha*max(V_rec,0))/tau



%while error(i) >= q_treshold
%end

end

