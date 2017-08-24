# Description
A recurrent model of V1 orientation tuning (Teich & Quian 2003) outfitted with a dynamic anti-Hebbian learning mechanism as described in:
Lange, Senden, Radermacher, and De Weerd (submitted). Learning multiple skills reveals competition rather than consolidation.

## Abstract
Well-trained skills are considered robust and long-lasting. Conforming to this, interference between different skills has only been reported when trained in single sessions separated by hours to days. Here, human participants underwent prolonged training in precisely judging two different orientations, chosen such that the two skills would rely on strongly overlapping neural population in early visual cortex, but would require different network connectivity. A network model of V1 orientation tuning fitted with a learning mechanism accurately predicted observed learning rates and skill levels, but also surprisingly strong interference among these highly trained skills during acquisition or retention. Prediction accuracy of the interference was maximized when no strict consolidation was assumed, such that all connections representing one skill could be modified by training another. Our data indicate a competition among skills for representation during which the skill with the longest training history wins at the cost of other, less-trained, skills.

## Files
This repository contains two files.
1. RM.m a MATLAB class implementation of the model.
2. LTI_Experiments.m: a script implementing the experiments reported in the aforementioned study.

This code is hosted at https://github.com/MSenden/LTI
The latest version may always be found here.

This software requires MATLAB R2015a or later. It was developed with access to the full suite of MATLAB add-on packages.
Some of these packages may be required to run the software.

# References
Teich AF, Qian N. Comparison Among Some Models of Orientation Selectivity. Journal of Neurophysiology. 2003; 89(4).

