---
# CastKin: An automated method for tracking head casting kinematics in freely moving _C. elegans_
CastKin is designed to extract the angular velocity of head casting in freely moving _C. elegans_. This package is designed to work with head centerlines obtained by [DeepLabCut](https://deeplabcut.github.io/DeepLabCut/README.html), but can be adapted for use with centerlines obtained from other popular methods (e.g., skeletonization, other pose trackers) as well! To get started, continue reading below.

## Features
- Jupyter Notebook-based for ease-of-use
- Utilizes [DeepLabCut](https://deeplabcut.github.io/DeepLabCut/README.html), a popular pose-tracking software, for extracting robust centerlines of the head
- Example DeepLabCut model weights, body part definitions, and skeleton are provided, but we _highly recommend_ training DeepLabCut on your own data for the best results.

## Contents
1. [Features](#features)
2. [Contents](#contents)
3. [Requirements](#requirements)
    - [Installation](#installation)
4. [Getting Started](#getting-started)
5. [Documentation](#documentation)
6. [Data Availability](#data-availability)
7. [References](#references)

## Requirements
We provide instructions for installing CastKin on Windows/MacOS/Linux below. While a GPU is **highly recommended** to use [DeepLabCut](https://deeplabcut.github.io/DeepLabCut/README.html), it is not required for use with CastKin.
You will need to download [Miniconda](https://www.anaconda.com/download/) and [Git](https://git-scm.com/install/) to set up the Python environment for CastKin and use the provided example notebook. We provide detailed instructions for how to do this below.

## Installation
1. Download [Miniconda](https://www.anaconda.com/download/) for your OS if you do not already have it (or Anaconda) installed.
2. Download [Git](https://git-scm.com/install/).
3. Open the program **Anaconda Prompt** to open a new terminal.
4. Enter the following command into your terminal to download this repository:
```bash
git clone https://github.com/lu-lab/CastKin.git
```
5. Navigate to the CastKin folder:
```bash
cd CastKin
```
6. Create the Conda environment:
```bash
conda env create -f environment.yml
```
7. Activate the environment:
```bash
conda activate castkin
```
8. Open Jupyter Notebook to run the example notebook:
```bash
jupyter notebook
```

## Documentation
Full documentation for CastKin is available at [https://erinshappell.github.io/castkin-docs/](https://erinshappell.github.io/castkin-docs/)

## Data Availability
Example DeepLabCut models and all videos used for validation may be found on PLACEHOLDER.

## References
CastKin uses [DeepLabCut](https://deeplabcut.github.io/DeepLabCut/README.html), an algorithm developed by [Mathis et al](https://doi.org/10.1038/s41593-018-0209-y). 

CastKin was written by Erin Shappell for Lu Lab, 2026.
