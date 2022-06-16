# Backward Simulation based on Random-graph Algorithms
SimRA (**Sim**ulation based on **R**andom-graphs **A**lgorithms) is a framework for simulating generic and complex evolutionary scenarios of multiple populations with subdivision and admixture. 
The algorithm generates a network structure, called Ancestral Recombination Graph (ARG), that models coalescence and genetic exchange events (recombinations) with SNPs. They support epistatic interactions 
between at most three loci at a time with user defined parameters. 

It has been developed in Java under the Eclipse Framework

# Pre-requisites

- Java JDK (version 1.6+)
- Apache Commons Math library. The jar file commons-math3-3.5.jar must be downloaded and stored in the SimRA\_lib folder. You can dowload it from [here] (http://commons.apache.org/proper/commons-math/download_math.cgi).
- Eclipse (Recommended, not required)

# Get Started
```sh
git clone git://github.com/ComputationalGenomics/SimRA
```

or you may want directly import the project on Eclipse (please see this [wiki](https://github.com/OneBusAway/onebusaway/wiki/Importing-source-code-into-Eclipse)). 
The jar file EpiSimRA.jar is in the subfolder called "binary".


# Citation

Please cite the following article if you use SimRA in your research:

A.P. Carrieri, F. Utro, L. Parida. Sampling ARG of multiple populations under complex configurations of subdivision and admixture. Bioinformatics, 2015.

if you use EpiSimRA, please cite: 

A. Bose, F. Utro, D.E. Platt, L.Parida. Multiway epistasis simulations for ARG sampling of multiple populations. (Submitted)
# Executable

For convenience, we also provide a precompiled version of SimRA usable via command line that you can download from the [binary folder](https://github.com/ComputationalGenomics/SimRA/tree/master/binary). 

Here we provide a short description on how to use it. For a more completed information please refer to the SimRA_UserManual.pdf file located in the binary folder.

## How to run?

```sh
$ java -jar EpiSimRa.jar [N] [r] [mu] [g] [iter] [eflag] [m] [s] [iter]
                 

 EXAMPLE:
         java -jar EpiSimRa.jar -N 10000 -r 0 -mu 1 -g 25 -m 10 20 30 40 -eflag 1 -s 0.3 0.3 0.3 -iter 100

 Above example is with selection on three loci and 4 randomly sampled extant units;

 As epistatic flag (eflag) is SET in the above example you need to provide four epistatic interaction coefficients less than s;
```

### Required parameters

```sh
- N : integer representing the population size
- r : double representing the recombination rate in cM/Mb/gen
- mu : double representing the mutation rate in  mut/bp/gen x 10^(-8)
- g : integer representing the segment length in Kb
- iter : Number of iterations
- eflag : Epistatic interaction flag
- m : array of integers less than N/2 representing the extant sample size
- s : array of doubles representing selection coefficient at each loci
- iter : number of iterations of the whole experiment
```

The parameter `s` when set to 0 executes the neutral case when there is no selection is in effect. It supports up 
to three-way epistatic interactions now with user defined values provided in the `s` parameter as shown in the example above. 

## Example
To execute sSimRA.jar it is necessary to reach the directory where the sSimRA.jar is located. 
Please refer to the example provided above to execute sSimRA.jar.


## Output

sSimRA generates a verbose output with stats on the ARG in the stats folder. The main output of ARG height is in
the text files with the following structure: 

`N~_m~_m0.2_g~_s~_r~.txt` 

More information about the ARG is stored in the files ending with `*STATS.txt`,`*L.txt` and `*S.txt`

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
