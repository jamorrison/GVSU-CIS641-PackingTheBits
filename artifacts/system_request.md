# Request Creator

- Jacob Morrison

# Project Sponsor

Dr. Erik Fredericks

# Business Need

This project exists to improve the efficiency of quality control of DNA sequencing experiments by merging multiple steps
of the process into a single executable.

# Business Requirements

- Merge coverage determination and statistical calculations
- Use multithreading to improve efficiency
- Integrate process within existing toolkit

# Business Value

I expect this project will reduce the time it takes to run the full quality control script, thereby decreasing the time
it takes to complete the initial analysis. Further, this reduction in time will benefit clinical users, who often work
in environments with limited resources and tight timelines.

# Special Issues or Constraints

- Project will be as a standalone tool before being integrated into existing toolkit
- Potential issue: Merging data from individual threads
- Potential issue: Allowing user to process a subset of data
