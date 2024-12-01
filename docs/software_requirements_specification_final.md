# Overview

This document serves as the culmination of the GVSU CIS 641 full semester project. It includes the final software
requirement specification, a change management plan, and traceability between requirements and various sundry diagrams.

# Software Requirements

Functional and non-functional requirements are provided. Additionally, a brief table of definitions is provided for
quick reference.

## Functional Requirements
### Command Line Interface Requirements
|  ID  | Requirement |
|:----:|:------------|
| FR01 | The command line interface shall accept an optional integer for the number of CPU threads. |
| FR02 | The command line interface shall accept an optional integer for the window step size. |
| FR03 | The command line interface shall accept an optional string for a prefix to all output file names. |
| FR04 | The command line interface shall accept an optional file path for the genomic windows with the top ten (10) percent of GC content. |
| FR05 | The command line interface shall accept an optional file path for the genomic windows with the bottom ten (10) percent of GC content. |
| FR06 | The command line interface shall accept a required file path for the CpG genomic locations. |
| FR07 | The command line interface shall accept a required file path for the input alignment file. |
| FR08 | The command line interface shall set a default value for the number of CPU threads when the corresponding command line option is not provided. |
| FR09 | The command line interface shall set a default value for the window step size when the corresponding command line option is not provided. |
| FR10 | The command line interface shall not prepend anything to the output file names when the prefix command line option is not provided. |
| FR11 | The command line interface shall check the number of CPU threads are valid greater than zero (0). |
| FR12 | The command line interface shall set the number of CPU threads to the default value when the given value is less than or equal to zero (0). |
| FR13 | The command line interface shall check the window step size is greater than zero (0). |
| FR14 | The command line interface shall set the window step size to the default value when the given value is less than or equal to zero (0). |

### File Input Requirements
|  ID  | Requirement |
|:----:|:------------|
| FR15 | The program shall read Sequence Alignment Map (SAM) files. |
| FR16 | The program shall read Binary Alignment Map (BAM) files. |
| FR17 | The program shall accept one SAM or BAM file as the alignment file input. |
| FR18 | The program shall extract chromosome sizes from the input SAM or BAM files. |
| FR19 | The program shall read Browser Extensible Data (BED) files. |
| FR20 | The program shall require BED files to be block gzipped. |

### File Output Requirements
|  ID  | Requirement |
|:----:|:------------|
| FR21    | The program shall have two output file types: a coverage distribution file (referred to as "covdist") and a coefficient of variation file (referred to as "cv"). |
| FR22    | The cv file shall have a two line header. |
| FR22.1  | The first line of the cv file shall say: "BISCUITqc Uniformity Table". |
| FR22.2  | The second line of the cv file shall say: "group mu sigma cv". |
| FR22.3  | The words on the second line of the cv file shall be separated by tab characters. |
| FR23    | All lines after the two line header in the cv file shall be tab separated. |
| FR24    | All lines after the two line header in the cv file shall contain (in the following order) a group name (see Functional Requirements FR24.1-FR24.12 for group names), the mean coverage, the coverage standard deviation, and the coverage coefficient of variation. |
| FR24.1  | The group name for coverage across bases covered by all applicable reads shall be "all_base". |
| FR24.2  | The group name for coverage across bases covered by reads with a mapping quality score (MAPQ) greater than or equal to forty (40) shall be "q40_base". |
| FR24.3  | The group name for coverage across CpGs covered by all applicable reads shall be "all_cpg". |
| FR24.4  | The group name for coverage across CpGs covered by reads with MAPQ greater than or equal to forty (40) shall be "q40_cpg". |
| FR24.5  | The group name for coverage across bases covered by all applicable reads and fall in the top GC windows shall be "all_base_topgc" when a top GC file is provided. |
| FR24.6  | The group name for coverage across bases covered by reads with MAPQ greater than or equal to forty (40) and fall in the top GC windows shall be "q40_base_topgc" when a top GC file is provided. |
| FR24.7  | The group name for coverage across CpGs covered by all applicable reads and fall in the top GC windows shall be "all_cpg_topgc" when a top GC file is provided. |
| FR24.8  | The group name for coverage across CpGs covered by reads with MAPQ greater than or equal to forty (40) and fall in the top GC windows shall be "q40_cpg_topgc" when a top GC file is provided. |
| FR24.9  | The group name for coverage across bases covered by all applicable reads and fall in the bottom GC windows shall be "all_base_botgc" when a bottom GC file is provided. |
| FR24.10 | The group name for coverage across bases covered by reads with MAPQ greater than or equal to forty (40) and fall in the bottom GC windows shall be "q40_base_botgc" when a bottom GC file is provided. |
| FR24.11 | The group name for coverage across CpGs covered by all applicable reads and fall in the bottom GC windows shall be "all_cpg_botgc" when a bottom GC file is provided. |
| FR24.12 | The group name for coverage across CpGs covered by reads with MAPQ greater than or equal to forty (40) and fall in the bottom GC windows shall be "q40_cpg_botgc" when a bottom GC file is provided. |
| FR25    | There shall be one covdist file for each set of coverages found (listed in Functional Requirements FR25.1-FR25.12). |
| FR25.1  | One covdist file shall be for bases covered by all applicable reads and shall have the tag "All Bases". |
| FR25.2  | One covdist file shall be for bases covered by all applicable reads with MAPQ greater than or equal to forty (40) and shall have the tag "Q40 Bases". |
| FR25.3  | One covdist file shall be for CpGs covered by all applicable reads and shall have the tag "All CpGs". |
| FR25.4  | One covdist file shall be for CpGs covered by all applicable reads with MAPQ greater than or equal to forty (40) and shall have the tag "Q40 CpGs". |
| FR25.5  | One covdist file shall be for bases covered by all applicable reads and fall in the top GC windows and shall have the tag "All Top GC Bases" when a top GC file is provided. |
| FR25.6  | One covdist file shall be for bases covered by all applicable reads with MAPQ greater than or equal to forty (40) and fall in the top GC windows and shall have the tag "Q40 Top GC Bases" when a top GC file is provided. |
| FR25.7  | One covdist file shall be for CpGs covered by all applicable reads and fall in the top GC windows and shall have the tag "All Top GC CpGs" when a top GC file is provided. |
| FR25.8  | One covdist file shall be for CpGs covered by all applicable reads with MAPQ greater than or equal to forty (40) and fall in the top GC windows and shall have the tag "Q40 Top GC CpGs" when a top GC file is provided. |
| FR25.9  | One covdist file shall be for bases covered by all applicable reads and fall in the bottom GC windows and shall have the tag "All Bot GC Bases" when a bottom GC file is provided. |
| FR25.10 | One covdist file shall be for bases covered by all applicable reads with MAPQ greater than or equal to forty (40) and fall in the bottom GC windows and shall have the tag "Q40 Bot GC Bases" when a bottom GC file is provided. |
| FR25.11 | One covdist file shall be for CpGs covered by all applicable reads and fall in the bottom GC windows and shall have the tag "All Bot GC CpGs" when a bottom GC file is provided. |
| FR25.12 | One covdist file shall be for CpGs covered by all applicable reads with MAPQ greater than or equal to forty (40) and fall in the bottom GC windows and shall have the tag "Q40 Bot GC CpGs" when a bottom GC file is provided. |
| FR26    | Each covdist file shall have a two line header. |
| FR26.1 | The first line of each covdist file shall say: "BISCUITqc Depth Distribution - TAG" where TAG is replaced by a tag unique to each covdist file (see Functional Requirements FR25.1-FR25.12 for the tags). |
| FR26.2 | The second line of each covdist file shall say: "depth count". |
| FR26.3 | The words on the second line of each covdist file shall be separated by tab characters. |
| FR27 | All lines after the two line header in the covdist file shall be tab separated. |
| FR28 | All lines after the two line header in the covdist file shall contain the coverage level in the first position of the line and the number of bases with that coverage level in the second position of the line. |

### Parallel Processing Requirements
|  ID  | Requirement |
|:----:|:------------|
| FR29   | The program shall use POSIX threads for multithreaded processing. |
| FR30   | The program shall split the chromosome into equal sized bins. |
| FR30.1 | The program shall allow the last bin to be less than the bin size if the chromosome cannot be broken up evenly. |
| FR31   | A thread shall store its block identification number. |
| FR32   | A thread shall store its bin start position. |
| FR33   | A thread shall store its bin end position. |
| FR34   | A thread shall store a bit array of positions in the bin that are covered by CpGs. |
| FR35   | A thread shall store a bit array of positions in the bin that are covered by top GC windows. |
| FR35.1 | The bit array shall be all zeroes (0s) if no top GC file is provided. |
| FR36   | A thread shall store a bit array of positions in the bin that are covered by bottom GC windows. |
| FR36.1 | The bit array shall be all zeroes (0s) if no bottom GC file is provided. |

### SAM/BAM File Processing Requirements
|  ID  | Requirement |
|:----:|:------------|
| FR37 | The program shall ignore unmapped reads. |
| FR38 | The program shall ignore reads that are not primary alignments. |
| FR39 | The program shall ignore reads that fail platform/vendor quality checks. |
| FR40 | The program shall ignore reads that are PCR or optical duplicates. |
| FR41 | The program shall ignore reads that are supplementary alignments. |
| FR42 | The program shall extract the CIGAR string from a mapped read alignment. |
| FR43 | The program shall increment one base along the reference if the CIGAR operation is an alignment match. |
| FR44 | The program shall not change the reference position if the CIGAR operation is an insertion. |
| FR45 | The program shall increment one base along the reference if the CIGAR operation is a deletion. |
| FR46 | The program shall increment one base along the reference if the CIGAR operation is a skipped region from the reference. |
| FR47 | The program shall not change the reference position if the CIGAR operation is a soft clipped base. |
| FR48 | The program shall not change the reference position if the CIGAR operation is a hard clipped base. |
| FR49 | The program shall increment one base along the reference if the CIGAR operation is a sequence match. |
| FR50 | The program shall increment one base along the reference if the CIGAR operation is a sequence mismatch. |
| FR51 | The program shall count a base as covered if the CIGAR operation is an alignment match, a sequence match, or a sequence mismatch. |
| FR52 | The program shall not count a base as covered if the CIGAR operation is an insertion, a deletion, a soft clipped base, a hard clipped base, or a skipped region from the reference. |

## Non-Functional Requirements
### Operating Requirements
|  ID   | Requirement |
|:-----:|:------------|
| NFR01 | The program shall work on Linux systems. |
| NFR02 | The program shall work on macOS systems. |
| NFR03 | The program shall run from a terminal. |
| NFR04 | The program shall work through a command line interface. |
| NFR05 | The program shall work on x86 CPU architectures. |
| NFR06 | The program shall work on ARM CPU architectures. |

### Performance Requirements
|  ID   | Requirement |
|:-----:|:------------|
| NFR07 | The program shall read a genomic window from a block gzipped BED file in two (2) seconds or less. |
| NFR08 | The program shall process at least one thousand (1000) reads per second. |
| NFR09 | The program shall process at least one million (1,000,000) basepairs of data per second. |
| NFR10 | The program shall use eight (8) gigabytes or less of random access memory. |
| NFR11 | The program shall compile in fifteen (15) seconds or less. |

### Usability Requirements
|  ID   | Requirement |
|:-----:|:------------|
| NFR12 | The program shall use English comments. |
| NFR13 | The program shall provide a help message for the user. |
| NFR14 | The program shall provide error messages upon system failure. |
| NFR15 | The program shall have a companion website describing usage. |
| NFR16 | The program shall provide an example dataset. |
| NFR17 | The program shall provide a tutorial for users. |

### Availability Requirements
|  ID   | Requirement |
|:-----:|:------------|
| NFR18 | The program shall be available for no cost. |
| NFR19 | The program shall be available to the public. |
| NFR20 | The program shall be available on GitHub. |
| NFR21 | The program shall be available in pre-compiled binaries. |
| NFR22 | The program shall be available through bioconda. |

### Maintainability Requirements
|  ID   | Requirement |
|:-----:|:------------|
| NFR23 | The program shall use Git for version control. |
| NFR24 | The program shall use semantic versioning. |
| NFR25 | Issues shall be acknowledged within two (2) business days. |
| NFR26 | Pull requests shall be acknowledged within two (2) business days. |
| NFR27 | Minor issues (e.g., bugs) shall be resolved within five (5) business days. |
| NFR28 | Major issues (e.g., feature requests) shall be resolved within fifteen (15) business days. |

## Abbreviation Definitions

| Abbreviation  | Definition |
|:--------------|:-----------|
| **SAM**       | Sequence Alignment Map file                                     |
| **BAM**       | Binary Alignment Map file                                       |
| **BED**       | Browser Extensible Data file                                    |
| **CIGAR**     | Compact Idiosyncratic Gapped Alignment Report string            |
| **Top GC**    | Top ten (10) percent of windows in the genome for GC content    |
| **Bottom GC** | Bottom ten (10) percent of windows in the genome for GC content |
| **MAPQ**      | Mapping quality score                                           |


# Change management plan
<Description of what this section is>

# Traceability links
<Description of this section>

## Use Case Diagram Traceability
| Artifact ID | Artifact Name | Requirement ID |
| :-------------: | :----------: | :----------: |
| UseCase1 | Move Player | FR5 |
| … | … | … |

## Class Diagram Traceability
| Artifact Name | Requirement ID |
| :-------------: |:----------: |
| classPlayer | NFR3, FR5 |
| … | … | … |

## Activity Diagram Traceability
<In this case, it makes more sense (I think, feel free to disagree) to link
to the file and to those requirements impacted>
| Artifact ID | Artifact Name | Requirement ID |
| :-------------: | :----------: | :----------: |
| <filename> | Handle Player Input | FR1-5, NFR2 |
| … | … | … |

# Software Artifacts
<Describe the purpose of this section>
* [I am a link](to_some_file.pdf)
