# Overview 

This document servers as the software requirements specification for the GVSU CIS 641 full semester project. Both
functional and non-functional requirements are described.

# Abbreviation Definitions

- **SAM file**: Sequence Alignment Map file
- **BAM file**: Binary Alignment Map file
- **BED file**: Browser Extensible Data file
- **CIGAR string**: Compact Idiosyncratic Gapped Alignment Report string
- **Top GC**: Top ten (10) percent of windows in the genome for GC content.
- **Bottom GC**: Bottom ten (10) percent of windows in the genome for GC content.
- **MAPQ**: Mapping quality score

# Functional Requirements
1. Command Line Interface Requirements
    1. The command line interface shall accept an optional integer for the number of CPU threads.
    2. The command line interface shall accept an optional integer for the window step size.
    3. The command line interface shall accept an optional string for a prefix to all output file names.
    4. The command line interface shall accept an optional file path for the genomic windows with the top ten (10)
    percent of GC content.
    5. The command line interface shall accept an optional file path for the genomic windows with the bottom ten (10)
    percent of GC content.
    6. The command line interface shall accept a required file path for the CpG genomic locations.
    7. The command line interface shall accept a required file path for the input alignment file.
    8. The command line interface shall set a default value for the number of CPU threads when the corresponding command
    line option is not provided.
    9. The command line interface shall set a default value for the window step size when the corresponding command line
    option is not provided.
    10. The command line interface shall not prepend anything to the output file names when the prefix command line
    option is not provided.
    11. The command line interface shall check the number of CPU threads are valid greater than zero (0).
    12. The command line interface shall set the number of CPU threads to the default value when the given value is less
    than or equal to zero (0).
    13. The command line interface shall check the window step size is greater than zero (0).
    14. The command line interface shall set the window step size to the default value when the given value is less than
    or equal to zero (0).
2. File Input Requirements
    1. The program shall read Sequence Alignment Map (SAM) files.
    2. The program shall read Binary Alignment Map (BAM) files.
    3. The program shall accept one SAM or BAM file as the alignment file input.
    4. The program shall extract chromosome sizes from the input SAM or BAM files.
    5. The program shall read Browser Extensible Data (BED) files.
    6. The program shall require BED files to be block gzipped.
3. File Output Requirements
    1. The program shall have two output file types: a coverage distribution file (referred to as "covdist") and a
    coefficient of variation file (referred to as "cv").
    2. The cv file shall have a two line header.
        1. The first line of the cv file shall say: "BISCUITqc Uniformity Table".
        2. The second line of the cv file shall say: "group mu sigma cv".
        3. The words on the second line of the cv file shall be separated by tab characters.
    3. All lines after the two line header in the cv file shall be tab separated.
    4. All lines after the two line header in the cv file shall contain (in the following order) a group name (see
    Functional Requirements 3.4.1-3.4.12 for group names), the mean coverage, the coverage standard deviation, and the
    coverage coefficient of variation.
        1. The group name for coverage across bases covered by all applicable reads shall be "all_base".
        2. The group name for coverage across bases covered by reads with a mapping quality score (MAPQ) greater than or
        equal to forty (40) shall be "q40_base".
        3. The group name for coverage across CpGs covered by all applicable reads shall be "all_cpg".
        4. The group name for coverage across CpGs covered by reads with MAPQ greater than or equal to forty (40) shall
        be "q40_cpg".
        5. The group name for coverage across bases covered by all applicable reads and fall in the top GC windows shall
        be "all_base_topgc" when a top GC file is provided.
        6. The group name for coverage across bases covered by reads with MAPQ greater than or equal to forty (40) and
        fall in the top GC windows shall be "q40_base_topgc" when a top GC file is provided.
        7. The group name for coverage across CpGs covered by all applicable reads and fall in the top GC windows shall
        be "all_cpg_topgc" when a top GC file is provided.
        8. The group name for coverage across CpGs covered by reads with MAPQ greater than or equal to forty (40) and
        fall in the top GC windows shall be "q40_cpg_topgc" when a top GC file is provided.
        9. The group name for coverage across bases covered by all applicable reads and fall in the bottom GC windows
        shall be "all_base_botgc" when a bottom GC file is provided.
        10. The group name for coverage across bases covered by reads with MAPQ greater than or equal to forty (40) and
        fall in the bottom GC windows shall be "q40_base_botgc" when a bottom GC file is provided.
        11. The group name for coverage across CpGs covered by all applicable reads and fall in the bottom GC windows
        shall be "all_cpg_botgc" when a bottom GC file is provided.
        12. The group name for coverage across CpGs covered by reads with MAPQ greater than or equal to forty (40) and
        fall in the bottom GC windows shall be "q40_cpg_botgc" when a bottom GC file is provided.
    5. There shall be one covdist file for each set of coverages found (listed in Functional Requirements 3.5.1-3.5.12).
        1. One covdist file shall be for bases covered by all applicable reads and shall have the tag "All Bases".
        2. One covdist file shall be for bases covered by all applicable reads with MAPQ greater than or equal to
        forty (40) and shall have the tag "Q40 Bases".
        3. One covdist file shall be for CpGs covered by all applicable reads and shall have the tag "All CpGs".
        4. One covdist file shall be for CpGs covered by all applicable reads with MAPQ greater than or equal to
        forty (40) and shall have the tag "Q40 CpGs".
        5. One covdist file shall be for bases covered by all applicable reads and fall in the top GC windows and shall
        have the tag "All Top GC Bases" when a top GC file is provided.
        6. One covdist file shall be for bases covered by all applicable reads with MAPQ greater than or equal to
        forty (40) and fall in the top GC windows and shall have the tag "Q40 Top GC Bases" when a top GC file is
        provided.
        7. One covdist file shall be for CpGs covered by all applicable reads and fall in the top GC windows and shall
        have the tag "All Top GC CpGs" when a top GC file is provided.
        8. One covdist file shall be for CpGs covered by all applicable reads with MAPQ greater than or equal to
        forty (40) and fall in the top GC windows and shall have the tag "Q40 Top GC CpGs" when a top GC file is
        provided.
        9. One covdist file shall be for bases covered by all applicable reads and fall in the bottom GC windows and
        shall have the tag "All Bot GC Bases" when a bottom GC file is provided.
        10. One covdist file shall be for bases covered by all applicable reads with MAPQ greater than or equal to
        forty (40) and fall in the bottom GC windows and shall have the tag "Q40 Bot GC Bases" when a bottom GC file is
        provided.
        11. One covdist file shall be for CpGs covered by all applicable reads and fall in the bottom GC windows and
        shall have the tag "All Bot GC CpGs" when a bottom GC file is provided.
        12. One covdist file shall be for CpGs covered by all applicable reads with MAPQ greater than or equal to
        forty (40) and fall in the bottom GC windows and shall have the tag "Q40 Bot GC CpGs" when a bottom GC file is
        provided.
    6. Each covdist file shall have a two line header.
        1. The first line of each covdist file shall say: "BISCUITqc Depth Distribution - TAG" where TAG is replaced by
        a tag unique to each covdist file (see Functional Requirements 3.5.1-3.5.12 for the tags).
        2. The second line of each covdist file shall say: "depth count".
        3. The words on the second line of each covdist file shall be separated by tab characters.
    7. All lines after the two line header in the covdist file shall be tab separated.
    8. All lines after the two line header in the covdist file shall contain the coverage level in the first position of
    the line and the number of bases with that coverage level in the second position of the line.
4. Parallel Processing Requirements
    1. The program shall use POSIX threads for multithreaded processing.
    2. The program shall split the chromosome into equal sized bins.
        1. The program shall allow the last bin to be less than the bin size if the chromosome cannot be broken up
        evenly.
    3. A thread shall store its block identification number.
    4. A thread shall store its bin start position.
    5. A thread shall store its bin end position.
    6. A thread shall store a bit array of positions in the bin that are covered by CpGs.
    7. A thread shall store a bit array of positions in the bin that are covered by top GC windows.
        1. The bit array shall be all zeroes (0s) if no top GC file is provided.
    8. A thread shall store a bit array of positions in the bin that are covered by bottom GC windows.
        1. The bit array shall be all zeroes (0s) if no bottom GC file is provided.
5. SAM/BAM File Processing Requirements
    1. The program shall ignore unmapped reads.
    2. The program shall ignore reads that are not primary alignments.
    3. The program shall ignore reads that fail platform/vendor quality checks.
    4. The program shall ignore reads that are PCR or optical duplicates.
    5. The program shall ignore reads that are supplementary alignments.
    6. The program shall extract the CIGAR string from a mapped read alignment.
    7. The program shall increment one base along the reference if the CIGAR operation is an alignment match.
    8. The program shall not change the reference position if the CIGAR operation is an insertion.
    9. The program shall increment one base along the reference if the CIGAR operation is a deletion.
    10. The program shall increment one base along the reference if the CIGAR operation is a skipped region from the
    reference.
    11. The program shall not change the reference position if the CIGAR operation is a soft clipped base.
    12. The program shall not change the reference position if the CIGAR operation is a hard clipped base.
    13. The program shall increment one base along the reference if the CIGAR operation is a sequence match.
    14. The program shall increment one base along the reference if the CIGAR operation is a sequence mismatch.
    15. The program shall count a base as covered if the CIGAR operation is an alignment match, a sequence match, or a
    sequence mismatch.
    16. The program shall not count a base as covered if the CIGAR operation is an insertion, a deletion, a soft clipped
    base, a hard clipped base, or a skipped region from the reference.

# Non-Functional Requirements
1. Operating Requirements
    1. The program shall work on Linux systems.
    2. The program shall work on macOS systems.
    3. The program shall run from a terminal.
    4. The program shall work through a command line interface.
    5. The program shall work on x86 CPU architectures.
    6. The program shall work on ARM CPU architectures.
2. Performance Requirements
    1. The program shall read a genomic window from a block gzipped BED file in two (2) seconds or less.
    2. The program shall process at least one thousand (1000) reads per second.
    3. The program shall process at least one million (1,000,000) basepairs of data per second.
    4. The program shall use eight (8) gigabytes or less of random access memory.
    5. The program shall compile in fifteen (15) seconds or less.
3. Usability Requirements
    1. The program shall use English comments.
    2. The program shall provide a help message for the user.
    3. The program shall provide error messages upon system failure.
    4. The program shall have a companion website describing usage.
    5. The program shall provide an example dataset.
    6. The program shall provide a tutorial for users.
4. Availability Requirements
    1. The program shall be available for no cost.
    2. The program shall be available to the public.
    3. The program shall be available on GitHub.
    4. The program shall be available in pre-compiled binaries.
    5. The program shall be available through bioconda.
5. Maintainability Requirements
    1. The program shall use Git for version control.
    2. The program shall use semantic versioning.
    3. Issues shall be acknowledged within two (2) business days.
    4. Pull requests shall be acknowledged within two (2) business days.
    5. Minor issues (e.g., bugs) shall be resolved within five (5) business days.
    6. Major issues (e.g., feature requests) shall be resolved within fifteen (15) business days.
