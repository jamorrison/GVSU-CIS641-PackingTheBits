# Packing the Bits

One common quality control analysis for DNA sequencing is finding the "coverage" of the genome. This shows the person
analyzing the data what fraction of the genome is covered by at least one "read" (a short fragment of sequenced DNA) and
how many reads, on average, span each base in the genome. In [BISCUIT](https://github.com/huishenlab/biscuit), an
aligner for whole genome bisulfite sequencing data, this is done through of combination of
[bedtools](https://github.com/arq5x/bedtools2) and GNU awk that is slow and highly repetitive. This project improves the
existing implementation through multithreaded reading of the SAM/BAM file and creating all metric and output files in
one command line interface call. Eventually, this will replace the current implementation in BISCUIT's quality control
pipeline.

## Team Members and Roles

* [Jacob Morrison](https://github.com/jamorrison/CIS641-HW2-Morrison)
  * Analyst, Designer, Developer, and Tester

## Prerequisites

* Operating system
  * Linux (tested on Alma Linux 9.4), OR
  * macOS (tested on Sonoma 14.6.1)
* Compilation tools
  * CMake (minimum version 3.21)
  * C compiler (tested with Apple Clang 15.0.0.15000309 and GCC 11.4.1)
  * Make
* System dependencies
  * cURL
  * zlib
  * POSIX threads
  * git
  * libdeflate (only needed for Linux)

## Run Instructions

### Installation

1. Clone repository:
```
git clone git@github.com:jamorrison/GVSU-CIS641-PackingTheBits.git
```
2. Move into repository directory and create and move to `build` directory:
```
cd GVSU-CIS641-PackingTheBits
mkdir build
cd build
```
3. Run CMake to configure build:
```
# Base configuration (release-mode build and default install path)
cmake ../

# Build in debug mode
cmake -DCMAKE_BUILD_TYPE=Debug ../

# Change installation path (replace <install-path> to desired location)
# Can also be combined with building in debug mode
cmake cmake -DCMAKE_INSTALL_PREFIX=<install-path> ../
```
4. Compile code
```
# Will not install binary anywhere
# (binary location is: build/src/coverage)
make

# Compile and install
# (binary location is: <install-path>/coverage)
make install
```

### Usage and Option Descriptions

```

Usage: coverage [options] <cpgs.bed.gz> <in.bam>

Options:
    -p STR    Prefix for output file names
    -b STR    Bottom 10 percent GC content windows BED file
    -t STR    Top 10 percent GC content windows BED file
    -s INT    Step size of windows [100000]
    -@ INT    Number of threads [1]
    -h        Print usage

```

* Required Arguments
  * __cpgs.bed.gz__: Path to block gzipped BED file of CpG locations in genome of interest
  * __in.bam__: Path to SAM or BAM file with aligned reads
* Optional Arguments
  * __-p__: Prefix (e.g., the sample name) for output files
  * __-b__: Path to block gzipped BED file of windows in the genome of interest with the bottom 10% of GC content
  * __-t__: Path to block gzipped BED file of windows in the genome of interest with the top 10% of GC content
  * __-s__: Window width for multithreaded processing (i.e., number of base pairs covered by each thread)
  * __-@__: Number of threads for processing SAM/BAM
  * __-h__: Print usage message

### Tutorial

Examples of using `coverage` can be found below. Example block gzipped BED files can be downloaded from GitHub
([link directly to file](https://github.com/jamorrison/GVSU-CIS641-PackingTheBits/releases/download/v0.1.0/bed_files.tar.gz))
or via the command line:
```
wget https://github.com/jamorrison/GVSU-CIS641-PackingTheBits/releases/download/v0.1.0/bed_files.tar.gz
```
The example BAM file is available in the `tests/` directory. All examples will be written as if the BED files and BAM
file are in the same directory, so please adjust paths as needed. Additionally, it will be assumed that `coverage` is in
your `PATH`.

#### Example 1: Basic usage

Find the coverage mean, standard deviation, and coefficient of variation (plus distribution of coverage values) for all
bases covered by all eligible reads, all bases covered by eligible reads with mapping quality score greater than or
equal to 40, plus CpGs covered by each of the two sets of reads used for all bases. This will use the default output
file names and only use one thread.
```
coverage cpgs.bgzip.bed.gz example.bam
```

#### Example 2: Multithreaded processing

To run with four threads, run:
```
coverage -@ 4 cpgs.bgzip.bed.gz example.bam
```
Replace `4` with your specific number of threads if you want to run with a different number of threads.

#### Example 3: Adding a prefix to output files

You can add an optional prefix to each of the output files:
```
coverage -p sample1 cpgs.bgzip.bed.gz example.bam
```
This would add `sample1_` (note the `_` is added by the program itself and does not need to be specified) to each of the
files output from `coverage`.

#### Example 4: Changing the window step size

By adjusting the window step size, you can balance the amount of memory used by each thread with the overall run time of
the tool. A larger step size will use more memory, but (in general) take less overall time due to fewer windows to
process. Conversely, a smaller step size will use less memory per thread, but take more time overall.
```
# Step size of one million (default is 100,000)
coverage -s 1000000 cpgs.bgzip.bed.gz example.bam

# Step size of one thousand
coverage -s 1000 cpgs.bgzip.bed.gz example.bam
```

#### Example 5: GC content regions

The `-t` and `-b` options will intersect the four default outputs with the BED files provided in those arguments and
create an additional four (or eight if both options provided) output files.
```
# Top 10% GC content windows only
coverage -t windows100bp.gc_content.top10p.bed.gz cpgs.bgzip.bed.gz example.bam

# Bottom 10% GC content windows only
coverage -b windows100bp.gc_content.bot10p.bed.gz cpgs.bgzip.bed.gz example.bam

# Both
coverage \
    -t windows100bp.gc_content.top10p.bed.gz \
    -b windows100bp.gc_content.bot10p.bed.gz \
    cpgs.bgzip.bed.gz \
    example.bam
```

#### Example 6: Putting it all together

Below is an example of all the options together:
```
coverage \
    -t windows100bp.gc_content.top10p.bed.gz \
    -b windows100bp.gc_content.bot10p.bed.gz \
    -@ 4 \
    -p sample1 \
    -s 1000000 \
    cpgs.bgzip.bed.gz \
    example.bam
```
