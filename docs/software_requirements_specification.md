# Overview 

This document servers as the software requirements specification for the GVSU CIS 641 full semester project. Both
functional and non-functional requirements are described.

# Abbreviation Definitions

- **SAM file**: Sequence Alignment Map file
- **BAM file**: Binary Alignment Map file
- **BED file**: Browser Extensible Data file
- **CIGAR string**: Compact Idiosyncratic Gapped Alignment Report string

# Functional Requirements
1. Command Line Interface Requirements
    1. The command line interface shall accept an integer for the amount of CPU threads.
    2. The command line interface shall accept a file path for the input.
    3. The command line interface shall check the inputs are valid.
    4. The command line interface shall accept a file path for the output.
    5. The command line interface shall accept an optional file of regions of interest.
    6. The command line interface shall accept an integer for chromosome bin width.
2. File Input Requirements
    1. The program shall read Sequence Alignment Map (SAM) files.
    2. The program shall read Binary Alignment Map (BAM) files.
    3. The program shall read Browser Extensible Data (BED) files.
    4. The program shall accept one SAM or BAM file as input.
    5. The program shall extract chromosome sizes from the input SAM or BAM files.
3. File Output Requirements
    1. The program shall output a BED-compliant file.
    2. The output shall be sorted first by chromsome, then by start position, then by end position.
    3. The output shall have a fourth (4th) column with the coverage value at that base.
    4. The output shall merge consecutive bases with the same coverage
    5. The output shall be tab-separated.
4. Parallel Processing Requirements
    1. The program shall use POSIX threads for multithreaded processing.
    2. The program shall split the chromosome into equal sized bins.
        1. The program shall allow the last bin to be less than the bin size if the chromosome cannot be broken up evenly.
    3. The program shall share an open read-only file handle of the input SAM or BAM file across threads.
    4. A thread shall store its bin start position.
    5. A thread shall store its bin end position.
    6. A thread shall store its output as a string.
5. SAM/BAM File Processing Requirements
    1. The program shall ignore unmapped reads.
    2. The program shall ignore reads that are not primary alignments.
    3. The program shall ignore reads that fail platform/vendor quality checks.
    4. The program shall ignore reads that are PCR or optical duplicates.
    5. The program shall ignore reads that are supplementary alignments.
    6. The program shall extract the CIGAR string from a mapped read alignment.
    7. The program shall increment one base along the read if the CIGAR operation is an alignment match.
    8. The program shall increment one base along the reference if the CIGAR operation is an alignment match.
    9. The program shall increment one base along the read if the CIGAR operation is an insertion.
    10. The program shall not change the reference position if the CIGAR operation is an insertion.
    11. The program shall not change the read position if the CIGAR operation is a deletion.
    12. The program shall increment one base along the reference if the CIGAR operation is a deletion.
    13. The program shall not change the read position if the CIGAR operation is a skipped region from the reference.
    14. The program shall increment one base along the reference if the CIGAR operation is a skipped region from the
    reference.
    15. The program shall increment one base along the read if the CIGAR operation is a soft clipped base.
    16. The program shall not change the reference position if the CIGAR operation is a soft clipped base.
    17. The program shall not change the read position if the CIGAR operation is a hard clipped base.
    18. The program shall not change the reference position if the CIGAR operation is a hard clipped base.
    19. The program shall increment one base along the read if the CIGAR operation is a sequence match.
    20. The program shall increment one base along the reference if the CIGAR operation is a sequence match.
    21. The program shall increment one base along the read if the CIGAR operation is a sequence mismatch.
    22. The program shall increment one base along the reference if the CIGAR operation is a sequence mismatch.
    23. The program shall count as base as covered if the CIGAR operation is an alignment match, a sequence match, or a
    sequence mismatch.
    23. The program shall not count as base as covered if the CIGAR operation is an insertion, a deletion, a soft clipped
    base, a hard clipped base, or a skipped region from the reference

# Non-Functional Requirements
1. Operating Requirements
    1. The program shall work on Linux systems.
    2. The program shall work on macOS systems.
    3. The program shall run from a terminal.
    4. The program shall work through a command line interface.
    5. The program shall work on x86 CPU architectures.
    6. The program shall work on ARM CPU architectures.
2. Performance Requirements
    1. The program shall read a region file in five (5) seconds.
    2. The program shall process at least one thousand (1000) reads per second.
    3. The program shall process at least one million (1,000,000) basepairs of data per second.
    4. The program shall use eight (8) gigabytes or less of random access memory.
    5. The program shall compile in fifteen (15) seconds or less.
3. Usability Requirements
    1. The program shall use English comments.
    2. The program shall provide a help message for the user.
    3. The program shall provide error messages upon system failure.
    4. The program shall having a companion website describing usage.
    5. The program shall provide an example dataset.
    6. The program shall provide a tutorial for users.
4. Availiability Requirements
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
