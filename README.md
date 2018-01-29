# Description
A recurrent model of V1 orientation tuning (Teich & Quian 2003) outfitted with a dynamic anti-Hebbian learning mechanism as described in:
Lange, Senden, Radermacher, and De Weerd (submitted). Interference with highly trained skills reveals  competition rather than consolidation.

## Abstract
Well-trained skills are considered robust and long-lasting. In line with this, interference between different skills has only been reported when trained in single sessions separated by hours to days. However, it is unknown whether, and to what extent, highly trained skills would resist interference. Here, human participants underwent prolonged training in precisely judging two different orientations, chosen such that the two visual skills would rely on strongly overlapping neural populations in early visual cortex, but would require different network connectivity. We hypothesized that through extensive training, a skill would be represented in an increasingly sparse manner, which could be achieved through interaction with read-out mechanisms selecting the most informative neurons for skill representation. This sparse representation would leave room for other skills to be represented within the same neural population, thus avoiding interference. Alternatively, even the most long-standing skills representing the highest expertise might be erased by other experience. We fitted a network model of V1 orientation tuning with learning rules according to these two hypotheses. We found that a visual skill trained extensively over a month could be erased by subsequent extensive training on a different visual skill. The prediction accuracy of the amount of interference was maximized when no consolidation was assumed, such that all connections representing one skill could be modified by training another. Thus, skills are not endowed by robustness over time thanks to consolidation. Instead, our empirical and model data indicate a competition among skills for representation during which the skill with the longest training history wins at the cost of other, less-trained, skills.

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

