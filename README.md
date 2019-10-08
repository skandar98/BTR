# Description
A recurrent model of V1 orientation tuning (Teich & Quian 2003) outfitted with a dynamic,
error-triggered, Hebbian learning mechanism as described in:
Lange, G., Senden, M., Radermacher, A., & De Weerd, P. (2020). Interfering with a memory without erasing its trace. Neural Networks, 121, 339–355. https://doi.org/10.1016/j.neunet.2019.09.027

## Abstract
Previous research has shown that performance of a novice skill can be easily interfered with by subsequent training of another skill. We address the open questions whether extensively trained skills show the same vulnerability to interference as novice skills and which memory mechanism regulates interference between expert skills. We developed a recurrent neural network model of V1 able to learn from feedback experienced over the course of a long-term orientation discrimination experiment. After first exposing the model to one discrimination task for 3480 consecutive trials, we assessed how its performance was affected by subsequent training in a second, similar task. Training the second task strongly interfered with the first (highly trained) discrimination skill. The magnitude of interference depended on the relative amounts of training devoted to the different tasks. We used these and other model outcomes as predictions for a perceptual learning experiment in which human participants underwent the same training protocol as our model. Specifically, over the course of three months participants underwent baseline training in one orientation discrimination task for 15 sessions before being trained for 15 sessions on a similar task and finally undergoing another 15 sessions of training on the first task (to assess interference). Across all conditions, the pattern of interference observed empirically closely matched model predictions. According to our model, behavioral interference can be explained by antagonistic changes in neuronal tuning induced by the two tasks. Remarkably, this did not stem from erasing connections due to earlier learning but rather from a reweighting of lateral inhibition.

## Data
The data related to the simulations is publicly available (https://doi.org/10.5061/dryad.6djh9w0wn).
Currently the data set is under review. In the meantine, the data can be dowloaded
from this temporary [link](https://datadryad.org/stash/share/RzgfKSZHyh_xaLqa6aEOVDHidDDyEpFgGi79_GlD1uU)

## Files
This repository contains two files.
1. RM.m a MATLAB class implementation of the model.
2. LTI_Experiments.m: a script implementing the experiments reported in the aforementioned study.

This code is hosted at https://github.com/ccnmaastricht/LTI
The latest version may always be found here.

This software requires MATLAB R2015a or later. It was developed with access to the full suite of MATLAB add-on packages.
Some of these packages may be required to run the software.

# References
Teich AF, Qian N. Comparison Among Some Models of Orientation Selectivity. Journal of Neurophysiology. 2003; 89(4).

