# Synergistic Antagonistic Interaction Detection
Implements the Synergistic Antagonistic Interaction Detection (SAID) framework.

# Abstract

There is abundant interest in assessing the joint effects of multiple exposures on human health.  This is often referred to as the mixtures problem in environmental epidemiology and toxicology.  Classically, studies have examined the adverse health effects of different chemicals one at a time, but there is concern that certain chemicals may act together to amplify each other's effects.  Such amplification is referred to as synergistic interaction, while chemicals that inhibit each other's effects have antagonistic interactions.  Current approaches for assessing the health effects of chemical mixtures do not explicitly consider synergy or antagonism in the modeling, instead focusing on either parametric or unconstrained nonparametric dose response surface modeling.  The parametric case can be too inflexible, while nonparametric methods face a curse of dimensionality that leads to overly wiggly and uninterpretable surface estimates. We propose a Bayesian approach that decomposes the response surface into additive main effects and pairwise interaction effects, and then detects synergistic and antagonistic interactions. A variable selection decision for each interaction component is also provided. This Synergistic Antagonistic Interaction Detection (SAID) framework is evaluated relative to existing approaches using simulation experiments and an application to data from NHANES.

ArXiv: https://arxiv.org/abs/2210.09279

# Install the R package (works on Linux, macOS, and Windows)

First install the 'devtools' package in R. If using Windows, install Rtools from here: https://cran.r-project.org/bin/windows/Rtools/rtools40.html before installing 'devtools'.

```
install.packages("devtools")
```
Next, install the package:

```
library(devtools)
install_github("shounakchattopadhyay/SAID")
```
