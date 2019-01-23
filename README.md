# Description
A recurrent model of V1 orientation tuning (Teich & Quian 2003) outfitted with a dynamic anti-Hebbian learning mechanism as described in:
Lange, Senden, Radermacher, and De Weerd (submitted). Interfering with a memory without disrupting its trace.

## Abstract
Synaptic consolidation is widely regarded to render well-trained skills long-lasting and resistant to disruption by subsequent learning. Support for this notion is mixed, however, as behavioral interference studies a) report inconsistent estimates of the time needed for memory traces to stabilize (from a few hours up to several days) and b) find interference to cease in some experimental conditions but persist in others. Here we integrate psychophysical experiments with computational modeling to definitively proof or reject the assumption that consolidation indeed protects skill memory traces from being overwritten. Human participants underwent extensive training in precisely judging two different orientations, chosen such that the two skills would rely on strongly overlapping neural populations in early visual cortex, but would require incompatible network connectivity. A network model of V1 orientation tuning fitted with a learning mechanism accurately predicted observed learning rates and skill levels, but also surprisingly strong interference among these highly trained skills during acquisition or retention. Prediction accuracy of the model was maximized when memory traces were assumed to remain labible, such that all connections representing one skill could always be modified by training on another skill. Our data indicate continuous competition among skills for representation such that the skill with the longest training history wins at the cost of other, less-trained, skills.

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

