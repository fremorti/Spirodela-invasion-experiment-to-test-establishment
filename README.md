# Invasion experiments test polyploid establishment and niche differentiation in greater duckweed Spirodela polyrhiza  
## Data and code
- analysis_comp_exps.R contains the annotated code to all statistical analysis and figure generation
- comp2022.csv and comp2023.csv contain the measured tetraploid and diploid nuclei counts from flow cytometry from both experiments
- cal.csv contains the calibration data from samples with known individual tetraploid proportion and measured tetraploid nuclei proportion
- dry weight.csv and dry weight2.csv contain the weekly dry weight measurements from both experiments
- envelopes.csv contains dry weights of two types of envelopes used in the experiment
- countweight contains dry weight measurements from samples with known number of individuals
- intrisic growth.csv contains the measured intrinsic growth rates that were separately measured and used to compare experimental results to in figure 2
- run_de_pop.sh is a shell script to call hpc_invasion_de_pop.R, which contains the ODE model fitting. The shell script can be used to run the computer intensive part on an external server.
