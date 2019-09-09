# Simulation based on Random-graph Algorithms
SimRA (**Sim**ulation based on **R**andom-graphs **A**lgorithms) is a framework for simulating generic and complex evolutionary scenarios of multiple populations with subdivision and admixture. 
The algorithm generates a network structure, called Ancestral Recombination Graph (ARG), that models coalescence and genetic exchange events (recombinations) with SNP.

An earlier version of our framework described in Carrieri et al (BioInformatics, 2015) and available here considered only the case of neutral scenario (where also STR polymorphisms can be generated). 
Recently, we extend the algorithm, sSimRA, that build multiple loci (_k_) selection, with multiway (_k_-way) epistasis for any arbitrary _k_, based on: (1) forward  (fwd-SimRA) and (2) backward (back-SimRA) simulation.


# Get Started
```sh
git clone git://github.com/ComputationalGenomics/SimRA
```

# Citation

Please cite the following article if you use sSimRA (backwards or forward) in your research:

A. Bose, F. Utro, D. Platt, L. Parida. sSimRA: multiple loci selection with multiway epistasis in coalescence with recombinations. (submitted)

and/or (for the neutral scenario): 

A.P. Carrieri, F. Utro, L. Parida. Sampling ARG of multiple populations under complex configurations of subdivision and admixture. Bioinformatics, 2015.


# Contributing

We welcome contributions but request that you follow the guidelines indicated [here](https://github.com/ComputationalGenomics/SimRA/blob/master/Contributing/Contributing.md).

# Apache License v. 2.0
Copyright 2015 IBM Corporation

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

# Contact

If you have suggestions, questions or comments regarding SimRA, please email us to: 

pa_ri_da (at) us.ibm.com  (remove the underscores)

# Open Source @ IBM

Find more open source projects on the [IBM Github Page](http://ibm.github.io/)
