# Backward Simulation based on Random-graph Algorithms
SimRA (**Sim**ulation based on **R**andom-graphs **A**lgorithms) is a framework for simulating generic and complex evolutionary scenarios of multiple populations with subdivision and admixture. The algorithm generates a network structure, called Ancestral Recombination Graph (ARG), that models coalescence and genetic exchange events (recombinations) with SNP as well as STR polymorphisms.

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

# Citation

Please cite the following article if you use SimRA in your research:

A.P. Carrieri, F. Utro, L. Parida. Sampling ARG of multiple populations under complex configurations of subdivision and admixture. Bioinformatics, 2015.

# Executable

For convenience, we also provide a precompiled version of SimRA usable via command line that you can download from the [binary folder](https://github.com/ComputationalGenomics/SimRA/tree/master/binary). 

Here we provide a short description on how to use it. For a more completed information please refer to the SimRA_UserManual.pdf file located in the binary folder.

## How to run?

```sh
$ java -jar SimRA.jar <whole path input directory> <name input file>.txt
<whole path output directory> <name output file> [-STR <STRs number> <state> <muSTR>]
```

### Required parameters

```sh
- <whole path input directory>: whole path of the directory when the input file is stored;
- <name input file>.txt: name for the input file that has to have a ".txt" extension;
- <whole path output directory>: whole path of the directory when the output files will be saved.
The output directory must be created by the user before executing SimRA;
- <name output file>: name for the output text files without extension;
```

### Optional parameters

```sh
- STR: it allows to get an output le containing information about STR polymorphisms;
- STRs number: number of STR loci - integer;
- state: initial state for each STR locus - integer;
- muSTR: STR mutation rate in mut/locus/gen x 10^(-4) - non negative real number;
```

## Example
To execute SimRA.jar it is necessary to reach the directory where the SimRA.jar is located. Some examples of command lines to execute SimRA.jar are the following:

### Without STR polymorphisms generation :
```sh
java -jar SimRA.jar ./scaffold_2admix.txt ./output scaffold_2admix
```
This command line allows to save the output files in the same directory where the executable SimRA.jar is located. In this case the input file is located in the same directory of SimRA.jar.

### With STR polymorphisms generation:

```sh
java -jar SimRA.jar ~/Desktop/inputSimRA/scaffold_2admix.txt ~/Desktop/outputSimRA/output
scaffold_2admix -STR 40 10 6.9
```
The whole path must be specified if the input/output directories inputSimRA, outputSimRA are located in a different directory of the executable SimRA.jar.

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
