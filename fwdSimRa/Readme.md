# Forward Simulation based on Random-graph Algorithms 
fwdSimRA (**Sim**ulation based on **R**andom-graphs **A**lgorithms) is a framework for simulating complex evolutionary scenarios of multiple populations
with recombination and coalescence with selection. The algorithm generates a network structure, called Ancestral Recombination Graph (ARG), that models
coalescence and recombinations with SNPs. This also takes into account the event when an individual's allele is under selection. 

It has been developed in C++ with the requirement of GSL and Openmp libraries. 

# Pre-requisites

 - C++ version 11 or later. 
 - GSL and Openmp libraries
 
# Get Started
```sh
git clone git://github.com/ComputationalGenomics/SimRA
```

# Executable 

You can recompile the fwdSimRA files using the following command 
```
g++ -std=c++11 fwdSimRA.cpp Individuals.cpp GlobalIndivs.cpp Events.cpp ChrInfo.cpp -o fwdSimRA -lgsl -fopenmp 
```
Some compilers might abort with an error about CBLAS linkage unavailability. In that case, use the extra -lgslcblas flag along with others. 

To run the package execute the following command: 
```
./fwdSimRA population_size(N) recombination_rate(r) mutation_rate(mu) chromosome_length(L) fitness_value(f) fitness_Value_of_first_allele (s1) fitness_value_of_second_allele (s2) Interaction_parameter (Delta) Number_of_interactions (iter)
```
For example if we want to run for N = 500 (250 males and 250 females) and r = 1e-7, mu = 15e-8, L = 20000, f = 0, s1 = 1, s2 = 0, delta = 0, iter = 100, we execute it as:
```
./fwdSimRA 250 0.0000001 0.00000015 20000 0 1 0 0 100
```
This will generate 4 files for m = {20, 40, 80, 120}

If you want to vary different m's, please edit the header of the fwdSSimRA.cpp file with the corresponding values of extant units(m). 

# Contributing

We welcome contributions but request that you follow the guidelines indicated [here](https://github.com/ComputationalGenomics/SimRA/blob/master/Contributing/Contributing.md).

# Apache License v. 2.0
Copyright 207 IBM Corporation

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
