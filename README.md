# AMC
Accurate mutation clustering from single-cell sequencing data.

## Requirements

* Linux systems.
* CMake3.0+.
* g++.
* [libeigen3](https://gitlab.com/libeigen/eigen)

## Installation

To build binary, do as follows:

```
tar -zxvf amc.tar.gz
cd amc
cmake .
make
```

The main program of AMC is generated in “bin” directory. Type following command if you want to add AMC to your system path:
```
make install
```

## Usage

AMC reasons mutation clusters from the genotype matrix derived from single-cell DNA sequencing.

Example:

```
amc -i testdata/example.txt -o testdata/example
```

## Input Files

### 1. Genotype Matrix

The SNVs of single cells are denoted as a genotype matrix D. Each row represents a mutation, and each column represents a single cell. Columns are separated by tabs. 

The genotype matrix can be ternary with element D[i,j] being

* 0 if mutation i is not observed in cell j,
* 1 if mutation i is observed in cell j, or
* 3 if the genotype information is missing

The genotype matrix can also be quaternary with element D[i,j] being

* 0 if mutation i is not observed in cell j,
* 1 if heterozygous mutation i is observed in cell j
* 2 if homozygous mutation i is observed in cell j
* 3 if the genotype information is missing

### 2. Cell names (optional)

A file listing the names of the single cells. Each row specifies the name of a single cell.
If no such file is specified, the cells are numbered from 1 to N (N is the number of cells).

### 3. Mutation names (optional)

A file listing the names of the mutations. Each row specifies the name of a mutation.
If no such file is specified, the mutations are numbered from 1 to M (M is the number of mutations).

## Output Files

### Recovered genotype matrix

The recovered genotype matrix is written to a file with suffix "recovered_gtm".

### Mutation assignments

The mutation assignments are written to a file with suffix "mutation_assignment".

### Input data for SCITE

The input file to run SCITE is created with suffix "scite_input".

### phylogenetic tree (optional)

When the path to the executable SCITE command is specified, the phylogenetic tree will be automatically reconstructed and written to a file in GraphViz format.

## Arguments

* `-i, --input <filename>` Replace \<filename\> with the file containing the genotype matrix.

* `-o, --output <string>` Replace \<string\> with the base name of the output file.

## Optional arguments

* `-c, --clabel <filename>` Replace \<filename\> with the path to the file containing the names of the cells.

* `-m, --mlabel <filename>` Replace \<filename\> with the path to the file containing the names of the mutations.

* `-S, --scite <String>` Replace \<String\> with the path to the executable SCITE command.

* `-K, --maxc <INT>` Set \<INT\> to a positive integer. This specifies the maximum number of mutation clusters to consider.

* `-a, --alpha <Double>` Set \<Double\> to fixed false positive rate.

* `-b, --beta <Double>` Set \<Double\> to fixed false negative rate.

* `-A, --max_alpha <Double>`  Set \<Double\> to the desired maximum false positive rate (effective when -a is not set). The default value is 0.05.

* `-B, --max_beta <Double>`  Set \<Double\> to the desired maximum false negative rate (effective when -b is not set). The default value is 0.5.

## Contact

If you have any questions, please contact zhyu@nxu.edu.cn.
