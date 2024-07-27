# RNN-AttractorNetwork-NewsomeShadlinTask

This repository contains MATLAB scripts to simulate recurrent neural networks (RNNs) that itinerate states given a Newsome and Shadlen style task. These scripts implement attractor-based neural firing models and explore various aspects of network dynamics, including bistability, multi-stable states, and the effects of different connectivity and noise parameters.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Scripts Overview](#scripts-overview)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)

## Installation

To run these scripts, you will need MATLAB installed on your machine. Ensure that you have the necessary toolboxes such as the Statistics and Machine Learning Toolbox.

## Usage

1. **Initialize and Set Parameters:** Modify the parameters within the script files as needed. Each script contains sections to set up the random number generator, parameters for single cell responses, connectivity parameters, and simulation parameters.

2. **Run the Scripts:** Execute the MATLAB scripts to simulate the neural network. For example, to run the attractor itinerancy simulation, use the following command in MATLAB:
   ```matlab
   run('ExampleCodes/attractor_itinerancy_graphic.m')
   run('Project/main.m')
   ```

3. **Visualize Results:** Each script generates graphical outputs to visualize the firing rates of excitatory cells and other relevant data. Ensure the graphics setup section is correctly configured to display or save the figures.

## Scripts Overview

### Main.m
The main script to run an attractor-based neural firing model on a stimulus regime designed to mimic the Newsome task. It includes pre-processing steps, general parameters, neural network setup, stimulus generation, and simulation execution.

### Subroutines
Contains helper functions to support the main scripts:
- `plotPermPopMeasures.m`
- `frProcess.m`
- `mapChoice.m`
- `populationAnalyses.m`

### ParameterExplorerScripts
Scripts to explore the parameter space for different aspects of the neural network:
- `PE_ExploreInputs.m`
- `PE_ExploreConnections.m`
- `PE_ExploreNoise.m`
- `PE_ExploreTauD.m`

### Classes
Defines classes to encapsulate neuron properties, input stimuli, and connection properties:
- `NeuronProperties.m`
- `InputStimulus.m`
- `ConnectionProperties.m`

## Project Structure

```
RNN-AttractorNetwork-NewsomeShadlinTask/
│
│
├── Project/
│   ├── Main.m
│   ├── Subroutines/ :: common subroutines
│   │   ├── plotPermPopMeasures.m
│   │   ├── frProcess.m
│   │   ├── mapChoice.m
│   │   └── populationAnalyses.m
│   │
│   ├── ParameterExplorerScripts/ :: simple grid search for parameters
│   │   ├── PE_ExploreInputs.m
│   │   ├── PE_ExploreConnections.m
│   │   ├── PE_ExploreNoise.m
│   │   └── PE_ExploreTauD.m
│   │
│   └── Classes/
│       ├── NeuronProperties.m
│       ├── InputStimulus.m
│       └── ConnectionProperties.m
│   
├── ExampleCodes/
│   ├── attractor_itinerancy_graphic.m
│   ├── longstim_rantau_spatial12.m
│   └── multiattractor_counting.m
│
└── README.md
```

