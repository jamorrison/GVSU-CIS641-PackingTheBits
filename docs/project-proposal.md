Team name: Packing the Bits

Team members: Jacob Morrison

# Introduction

(In 2-4 paragraphs, describe your project concept)

Outline
- Background
  - Genomics background
  - Current tools used
  - Other possibilities
  - Why I'm writing my own
- Project Description
  - Written in C to plug into BISCUIT
  - Multithreaded to speed up coverage calculations
  - Calculate output metrics on the fly
- Future Use
  - Plug into BISCUIT after class
  - Probably not enough for it's own tool, but could theoretically be published as a standalone tool

# Anticipated Technologies

- OS: macOS or Linux
- Editor: vim
- Language: C
- Third-Party Libraries:
  - [htslib](https://github.com/samtools/htslib)
  - POSIX threads
  - [zlib](https://www.zlib.net)

# Method/Approach

(What is your estimated "plan of attack" for developing this project)

# Estimated Timeline

(Figure out what your major milestones for this project will be, including how long you anticipate it *may* take to reach that point)

# Anticipated Problems

Potential problems include, but are not limited to:

- Analysis
- Design
  - Creating class diagrams
- Implementation
  - Handling reads that span adjacent regions handled by different threads
  - Merging regions of similar coverage that span adjacent regions
- Testing
  - Unit testing
