#SimRA Executable
SimRa is also provided as a Java executable via command line.


#How to run?

Here we provide a short description on how to use it. For a more completed information please refer to the SimRA_UserManual.pdf file located in this folder.


```sh
$ java -jar SimRa.jar <whole path input directory> <name input file>.txt <whole path output directory> <name output file> [-STR <STRs number> <state> <muSTR>]
```

##Required parameters

```sh
- <whole path input directory>: whole path of the directory when the input le is stored;
- <name input file>.txt: name for the input le that has to have a ".txt" extension;
- <whole path output directory>: whole path of the directory when the output les will be saved. The output directory must be created by the user before executing SimRA;
- <name output file>: name for the output text les without extension;
```

##Optional parameters

```sh
- STR: it allows to get an output le containing information about STR
polymorphisms;
- STRs number: number of STR loci - integer;
- state: initial state for each STR locus - integer;
- muSTR: STR mutation rate in mut/locus/gen x 10^􀀀4 - non negative real
number;
```

#Example
To execute SimRA.jar it is necessary to reach the directory where the SimRA.jar is located. Some examples of command lines to execute SimRA.jar are the following:

##Without STR polymorphisms generation:
```sh
java -jar SimRA.jar ./scaffold 2admix.txt ./output scaffold 2admix
```
This command line allows to save the output files in the same directory where the executable SimRA.jar is located. In this case the input file is located in the same directory of SimRA.jar.

##With STR polymorphisms generation:

```sh
java -jar SimRA.jar ~/Desktop/inputSimRA/scaffold 2admix.txt ~/Desktop/outputSimRA/output scaffold 2admix -STR 40 10 6.9
```
The whole path must be specified if the input/output directories inputSimRA, outputSimRA are located in a different directory of the executable SimRA.jar.

