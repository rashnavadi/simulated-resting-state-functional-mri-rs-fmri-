**Resting-State fMRI Dynamic Functional Connectivity Analysis**

**Overview**

This repository contains code and data for testing the performance of two dynamic functional connectivity analysis approaches on resting-state fMRI (rs-fMRI) data. Specifically, we compare:

Hierarchical Observational COnnectivity (HoCo)

Sliding Window Cross-Correlation (SWC)

Each method is combined with k-means clustering to assess their sensitivity in detecting state transitions in simulated rs-fMRI data.

**Simulation Details**

To evaluate the performance of these approaches, we generated simulated rs-fMRI data with predefined state transitions. 
The simulation involves three pure sine wave signals (x, y, and z) with the following characteristics:
x, y, z Signals: Simulated as pure sine waves with specifications matching typical rs-fMRI data, but with minimal noise added.
State Transitions: Imposed two states on the network with a transition occurring at the 150-second (or 75 timepoint with TR=2 sec) mark.

State 1 (0 to 150 seconds):
Correlation between x and y: 0.4
Correlation between x and z: 0.8
Correlation between y and z: 0.8

State 2 (150 to 300 seconds):
Correlation between x and y: 0.8
Correlation between x and z: 0.4
Correlation between y and z: 0.8

These simulations were repeated for 1000 iterations to generate robust data sets for analysis.

**No Noise Analysis: hard-coded timeseries of x, y and z**
The analysis is performed in the script no_noise.m. We expected both methods to detect the transition at the exact time of 150 seconds. Here, I used three sine timeseries and hard-coded that phase shifts between the x, y and z timeseries to impose the correlations.

The results confirmed this expectation:

HoCo + k-means detected the transition at timepoint 76.
SWC + k-means detected the transition at timepoint 51.

Considering the repetition time (TR) is 2 seconds, these detections correspond to the ground truth transition time of 150 seconds for both methods.

**generateCorrelatedSinusiods.m: a general version for simulating no noise x, y and z timeseris where I coded how to find the timeseries with desired correlaitons between the timeseries**

**in the script called call_generateCorrelatedSinusoids.m I iterated trhough the timesereis generation for 1000 timese and run hoco/swc + kmeans on each iteration to compare their performance**

still in this analyses we expect to have a perfect detection of state transitin from both methods, which we did.

**Enhanced Simulation with Noise: apply_diff_noise_levels.m**

In the next step of the analysis, we aimed to make the x, y, and z time series more similar to actual rs-fMRI time series by adding pink noise and Gaussian noise while maintaining the correlation conditions between them. The Signal-to-Noise Ratio (SNR) for rs-fMRI is defined as 50, which is defined as the mean of the time series over time divided by the standard deviation of the entire noise for each time series. To achieve this:

The mean of the time series was rescaled to 10000, while keeping the SNR at 50.
This generation of time series was repeated 1000 times.

The performance of HoCo + k-means and SWC + k-means was compared again to check if they can detect the imposed transition of states and which method has higher sensitivity.


**Repository Structure**
code/: Scripts for generating simulations, performing HoCo and SWC analyses, and applying k-means clustering. The x, y and z timeseries are generated in the code.

output/: Specific branch dedicated to storing analysis outputs.

Branches

main: Contains the primary codebase and data.

outputs: Dedicated branch for storing analysis outputs. This branch contains all the results generated by the scripts in the code/ directory.

**Contributions**

Contributions are welcome! Please create a new branch for any feature additions or bug fixes and submit a pull request for review.

**Contact**

For any questions or further information, please contact Tahereh (Tara) Rashnavadi@: tahereh.rashnavadi@gmail.com

