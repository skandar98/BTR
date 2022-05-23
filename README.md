
#BTR
Nonlinear and linearised recurrent models of V1 orientation tuning

This repository contains all resources used in the Bachelor Thesis: 

"Changes in System Dynamics Reflect Consolidation of Memory during Formation and Interference"

## Abstract:

Recent studies have proposed and tested multiple neural network models to capture
some characteristics of perceptual learning in the primary visual area, such as specificity, transfer and interference between visual memory traces. The present thesis aims
to complement these results by qualitatively describing the behaviour of the nonlinear
recurrent model used by Lange, Senden, Radermacher, and De Weerd (2020). We replicated their simulations of an interference condition, and assessed the dynamics of the
model throughout the network’s training by computing the kinetic energy, and applied
principal component analysis to identify the terminal states of the system after each
trial. We present evidence that two principal components (PC1 and PC2) can capture
the system’s end-states, where PC2 was associated with the system’s decision, and PC1
was associated with the kinetic energy of the system. Throughout the training, for each
reference angle we observed a decrease in system velocity, and a divergence between terminal states along PC1. Similar results were obtained from a linearised version of the
model, which however traded accuracy for computational efficiency. Based on our results, we suggest that anti-Hebbian learning in recurrent networks increases performance
by further separating equilibrium points in phase-space, thereby allowing the system to
reach a terminal state that is closer to the equilibrium point associated with the correct
decision

## Contents:

- In BTR/RNN/:
1. The Matlab classes RM_GN.m and RM_GL.m, the nonlinear and linear recurrent models.
2. The Matlab scripts main_nonlin.m and main_lin.m, used to run the simulations, collect and visualise data.
3. The MAtlab script simulate_forward.m used for confirmng the validity of the utilized simulation timespan.
4. old/ contains further recurrent model versions that were ommited from the final paper.

- BTR/FPAnalysis/ contains depricated files containing additional analyses that were ommited from the final paper. 
 
