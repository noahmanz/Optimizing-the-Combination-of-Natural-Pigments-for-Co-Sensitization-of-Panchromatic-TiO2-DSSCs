# Optimizing-the-Combination-of-Natural-Pigments-for-Co-Sensitization-of-Panchromatic-TiO2-DSSCs

  This repository contains the source code, raw data and generated figures used in the MS thesis “Optimizing the Combination of Natural Pigments for Co-Sensitization of Panchromatic TiO2 Dye Sensitized Solar Cells” by Noah Manz. The thesis is published online at:

https://www.proquest.com/openview/5e34518343751bc814c51ea0720afd66/1?pq-origsite=gscholar&cbl=18750&diss=y

  The included “Optimization_Script.py” file allows users to obtain optimized light harvesting efficiency (LHE) spectra for an arbitrary set of constituent dyes by generating analogs of the “Empirical_Dye_Solutions_Volume_Fractions.csv” and “UVVIS_Absorbance_Anode_Adsorbed” or “UVVIS_Absorbance_Bulk_Solution.csv” files.


# About the Repository:

  The premise of this work is to use empirical UV/VIS absorbance data to construct higher-dimensional surfaces which interpolate the measured data. This has been achieved with Radial Basis Function (RBF) Interpolation from the Scipy Interpolate package. For example, given a set of 6 constituent dyes, this would be a surface in R6 which interpolates all measured absorbance values at a given wavelength. Assuming that absorbance data is available for C combinations of these 6 constituent dyes, this would be a function interpolating C points in R6. This process is repeated for all wavelengths for which data is available. Then, each interpolation function can be sampled and concatenated to yield a spectrum. These interpolation functions can be sampled at any point (a point in R6) to produce estimates of the absorbance profile for that particular combination. As C increases, the accuracy of the model, and therefore the accuracy of this estimate, also increase. This allows for evaluation of potentially thousands of dye combinations from a more modest empirical database.

  To learn more about Radial Basis Function interpolation and why it was used in this work (in lieu of a simple rule-of-mixtures), check out the following video on my YouTube channel:

<a href="http://www.youtube.com/watch?feature=player_embedded&v=KSHNrELYn9g
" target="_blank"><img src="http://img.youtube.com/vi/KSHNrELYn9g/0.jpg" 
alt="Radial Basis Function Interpolation" width="240" height="180" border="10" /></a>

  To obtain absorbance-optimized dye solutions, the results of these interpolation functions are compared with the AM1.5G solar irradiance spectrum and evaluated in one of three so-called “fitment conditions”. Three fitment conditions (as well as RBF interpolation functionality) are hard-coded into the “Optimization_Script.py” file (i.e Pearson Correlation, Integral Value & Covariance Value). These are conditions which evaluate the commensurability of a particular LHE spectrum with AM1.5G, albeit with slightly different assumptions. This work has been performed hypothesizing that dye combinations maximizing one of these fitment conditions will also maximize DSSC performance.

  To learn more about the different fitment conditions used in this work, check out the following video:

<a href="http://www.youtube.com/watch?feature=player_embedded&v=D9Z7w32d_Ts&t
" target="_blank"><img src="http://img.youtube.com/vi/D9Z7w32d_Ts/0.jpg" 
alt="Fitment Conditions for AM1.5G Commensurability" width="240" height="180" border="10" /></a>

# How to Use this Code:

  To perform this optimization, download the .py file, UVVIS_Absorbance files, and “Empirical_Dye_Solutions_Volume_Fractions.csv” file. This code can be run “as-is” and will replicate the published results, however, proprietary data can also be substituted for optimization of new combinations, and the AM1.5G spectrum can also be subtitled for evaluation of DSSC performance in extra-terrestrial conditions. To do this, ensure that your data is formatted in the same manner as in the provided data files. In particular, that the UVVIS data is column-oriented (with each unique dye combination populating a column) where each row represents a particular wavelength. You may include as many combinations as you’d like (See example table below). Ensure however that the wavelength domain is constant for all spectra. This has been set to 300 - 800 nm (by 1 nm). To adjust this, see Line 36 in the .py file. 
  
| Wavelength (nm) | Combo. 1 Abs. | Combo. 2 Abs. | ... | Combo. C Abs. |
| --- | --- | --- | --- | --- |
| 300 | 0.1 | 0.2 | ... | 0.3 |
| 301 | 0.4 | 0.5 | ... | 0.6 |
| ... | ... | ... | ... | ... |
| 800 | 0.7 | 0.8 | ... | 0.9 |
  
  Also populate an analog of the “Empirical_Dye_Solutions_Volume_Fractions.csv” file. This file contains the volume fractions of constituent dyes (the combinations) which produced the data in the UVVIS file. This file is row-oriented, so ensure that the transpose of the header row in the UVVIS file produces the first column in the “Empirical_Dye_Solutions_Volume_Fractions.csv” file. In addition to evaluating any number of combinations, C, you may also populate the Volume Fraction file to include as many constituents as you’d like. If this is more than 6  constituents however, you will have to adjust some pieces of the code. Instructions for how to do this have been written into the script showing a fictitious evaluation of 7 constituent dyes. These rules generalized for D number of constituents. If you are evaluating less than 6 constituent dyes, it is easier to just set the volume fractions of any remaining dyes to 0.
  
| Dye Combination | V[f] of Dye 1 | Vf of Dye 2 | ... | Vf of Dye 6 |
| --- | --- | --- | --- | --- |
| Combo. 1 | 0.1 | 0.2 | ... | 0.3 |
| Combo. 2 | 0.4 | 0.5 | ... | 0.6 |
| ... | ... | ... | ... | ... |
| Combo. C | 0.7 | 0.8 | ... | 0.9 |
 

  By default, this code will output 6 different pieces of information. It will print the resolution of the volume fraction array and number of combinations being evaluated in the console. It will also print the combinations of dyes that maximize each of the 3 fitment conditions. It will output a plot showing the AM1.5G solar irradiance spectrum, and it will output 2 plots for each fitment condition. The latter of these plots shows the LHE spectrum of the optimized dye vs. AM1.5G, and the former shows the value of the fitment condition associated with each dye combination evaluated. You can prevent any of these plots from being generated by setting allow = False in each plot function.

# Folders Containing Plots:

  The plots contained in the “IV Plots” and “UV_VIS Plots” folders are intended as a convenience for those with limited plotting or Python skills. They may be reproduced with express permission only (which I am happy to grant). You can also email me directly at: nmanzf35@gmail.com
