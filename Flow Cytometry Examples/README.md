# Inverter characterization experiment

The data in this directory is a selected subset from an experiment conducted in May, 2014 by Swati Banerjee Carr as part of her Ph.D. thesis, providing preliminary data towards the following publications:

* Swati Banerjee Carr, "Reliable Gene Expression and Assembly for Synthetic Biological Devices in E. Coli Through Customized Promoter Insulator Elements and Automated DNA Assembly." Ph.D. Thesis, Boston University, July, 2016.
* Swati Banerjee Carr, Jacob Beal, and Douglas M. Densmore. "Reducing DNA context dependence in bacterial promoters." PLoS ONE 12(4): e0176013, April 2017

## Experiment Description

Circuit: Arabinose-induced TetR inverter RBS-variations

Selected Circuit:

- INV-155: pBAD-RBS30-GFP-Term::pConst(J23100)-RBS30-AraC-Term::pBAD-RBS30-TetR-Term::pTetR-RBS34m1-RFPm-Term

Process Controls:

- SpheroTech RCP-30-5A rainbow calibration particles
- Bioline competent cells (unmodified)
- Constitutive GFP
- Constitutive RFP
- Constitutive GFP + RFP

Growth Conditions:

- Biological triplicates of each inverter were grown overnight in 2X strengh LB Broth with Ampicillin. 

- Each inverter replicate was diluted (5uL in 145uL) at 6PM into LB broth containing Ampicillin and L-Arabinose (8 L-Ara concentration titration) 

- Inverters were allowed to grow until 8AM the following morning, when they were diluted into PBS (1/100X dilution) and run through the flow cytometer. 
 
## Example Analyses

Example analyses were prepared using Matlab and TASBE Flow Analytics 8.3.2.

The analyses in this directory should be executed in the following order

- `colormodel/make_color_model.m`: This script uses the process controls to set up conversion from raw flow cytometry data to properly compensated and calibrated data. The conversion is contained in the `ColorModel` object in `CM140521.mat`
- `analysis/batch_analysis.m`: This script processes the experimental data, converting the raw files, organizing by replicates, and outputting full statistics in `TetR-batch.mat` and extracts in the CSV directory.
- `plot_inverter.m`: This script plots the results of the experiment as a summary transfer curve, from which model parameters may be extracted.
