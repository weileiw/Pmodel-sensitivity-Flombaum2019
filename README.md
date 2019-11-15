# Pmodel-sensitivity-Flombaum2019

This repository contains the sensitivity test code used in Flombaum et al 2019 (Nature Geoscience -- Picophytoplankton lineages display clear niche partitioning but overall positive response to future ocean warming).

The main driver is "SENSI_TEST.m", which sets up the model environment and defines model parameters.

The main function then calls eqPcycle_sensi.m for calculation of new field of phosphorus (DIP, DOP, and POP) based on altered model parameters.


unitlity scripts:

nsnew.m C.T. Kelly's Newton Armijo solver

mfactor.m Timothy A. Davis' LINFACTOR VERSION 1.1.0, Nov 1, 2007 Copyright 2007, Timothy A. Davis, University of Florida

d0.m makes a sparse diagonal matrix given a 3d field
