Team name: Packing the Bits

Team members: Jacob Morrison

# Introduction

DNA sequencing is a common practice in biomedical research. The sequencing process requires collecting DNA from a
sample (e.g., a tumor biopsy or a healthy tissue nearby), extracting the DNA from the cells, and then fragmenting the
full length DNA strands into smaller chunks. These chunks are then placed on a "sequencer," where the individual bases
are read from the DNA fragments. The output from the sequencer is called a "read" and is then aligned back to a
reference genome (the human reference genome is about three billion bases longer) in order to find where the small
fragment (usually around 100-150 bases in length) came from.

One possible analysis that can be done is to find the "coverage" across the entire genome. To do this, you count the
total number of reads that covers each individual loci, or base, in the genome. One example of performing this analysis
is as a part of a quality control pipeline within [BISCUIT](https://github.com/huishenlab/biscuit) (a project I maintain
for work) to verify a sequencing experiment consistently covered the entire genome. Currently, the pipeline uses
[bedtools](https://github.com/arq5x/bedtools2) to do this, but it is a slow process. One option for improving the speed
of this process would be to use [mosdepth](https://github.com/brentp/mosdepth), a tool specifically designed for finding
coverages. However, it does not calculate the exact metrics that I need.

I propose writing a small tool that improves upon the existing bedtools implementation and calculates the metrics I need
on the fly. The tool will be written in C for plugging in to BISCUIT following the completion of class. It will also be
multithreaded to improve the efficiency of the single-threaded bedtools option. Depending on the final state of the
project, it may be sufficient to warrant an academic publication (either as a pre-print or a peer-reviewed application
note).

# Anticipated Technologies

- OS: macOS or Linux
- Editor: vim
- Language: C
- Third-Party Libraries:
  - [htslib](https://github.com/samtools/htslib)
  - POSIX threads
  - [zlib](https://www.zlib.net)

# Method/Approach

Broadly, I plan to start off with doing as much analysis and design ahead of time as I can. This will include defining
as many requirements at the outset of the project (although I anticipate having to revise these during the
implementation stage). I plan to give myself the most time on the implementation stage, with an approximately equal
amount of time for analysis, design, and implementation. Testing will have a shorter amount of dedicated time; however,
I plan to write some of the tests during the implementation stage. For the sake of this project, most of my tests will
focus on the output of the tool to ensure it returns the same values as before. I will include some unit tests on some
of the smaller functionality, but these won't be the focus of my tests.

# Estimated Timeline

| Task | Estimated Completion Date |
|:-----|:--------:|
| Analysis       | 30 September 2024 |
| Design         | 14 October 2024   |
| Requirements   | 14 October 2024   |
| Implementation | 14 November 2024  |
| Testing        | 31 November 2024  |
| Wrapping Up    | 6 December 2024   |

# Anticipated Problems

Potential problems include, but are not limited to:

- Design
  - Creating class diagrams
- Implementation
  - Handling reads that span adjacent regions handled by different threads
  - Merging regions of similar coverage that span adjacent regions
  - Handling regions of interest to calculate metrics for
- Testing
  - Unit testing
