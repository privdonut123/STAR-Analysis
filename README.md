# TSSA for pp&uarr; 510 GeV collsions at STAR

This project contains scripts and data for hadron analysis, test analysis, and other related tasks.

Much of the code is experimental and is not used in the final analysis. The final analysis macro is contained in `hadron_Analysis.C`.

## Directory Structure

- `hadron_Analysis.C`: Contains the main script for hadron analysis.
- `test_Analysis.C`: Contains the test script for the analysis.
- `HelloWorldExample/`: Contains example scripts and data.
- `Analysis.C`: Contains the main analysis script.
- `README.md`: This file.
- `Run_17_polariztion.txt`: Contains data for Run 17 polarization. This file isn't used in the analysis, but Run 22 polarization is used instead.
- `daq/`: Contains scripts for FTS (Foward Tracking System). 🚩
- `env.sh`: Environment setup script for FTS. Used for debugging. Not used in the analysis. 🚩
- `fitProjections.C`: Script for fitting ECAL projections onto HCAL.
- `include/`: Contains included header files for colored output.
- `input/`: Contains input data files.
- `output/`: Contains output data files.
- `polarization.C`: Script for polarization analysis. Currently not used.
- `read.sh`: Script for running analysis. Currently not used. 🚩
- `spins/`: Contains data for spins for Run 17. Currently not used. 🚩
- `temp_gccflags.c`: Temporary file for GCC flags.
- `time_stamp.C`: Script for timestamping. Currently not used. 🚩

## How to Run

1. Clone the repository.
2. From the terminal, run the following commands:

```
cd STAR-Analysis
stardev
cons
```
Running `cons` will build all of the neccessary libraries for the analysis. This may take a while. Once the libraries are built, run any of the analysis macros, depending on what you want to do.

For example:
```
root4star -b -q -l hadron_Analysis.C
```

All of the output files will be stored in the `output/` directory.