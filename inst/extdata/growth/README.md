# Data sets for GM and GDM analysis

## Source

Susan M. Motch Perrine, Tim Stecko, Thomas Neuberger, Ethylin W. Jabs, 
Timothy M. Ryan and Joan T. Richtsmeier, 2007. 
Integration of brain and skull in prenatal mouse models of 
Apert and Crouzon Syndromes. **Front. Hum. Neurosci.**, 
URL: https://doi.org/10.3389/fnhum.2017.00369

## Files

The data files are:

- `CZEM_mut_global.xyz`: Crouzon mutant embryonic mouse age E17.5
- `CZEM_wt_global.xyz`: Crouzon unaffected littermate embryonic mouse age E17.5
- `CZP0_mut_global.xyz`: Crouzon mutant newborn mouse 
- `CZP0_wt_global.xyz`: Crouzon unaffected littermate newborn mouse

## Landmarks and data subsets for analysis

The subset I used includes landmarks that are located all over the 
neurocranium (not the facial skeleton). From the 43 landmarks in these file, 
I chose 10 landmarks used in the attached publication:
 
- numbers: 1  2 13 17 19 20 33 38 40 41
- abbreviations: amsph   bas  loci  lpto  lsqu  lsyn  roci  rpto  rsqu  rsyn

amsph and bas are on the midline. The other landmarks are bilateral. 
Landmarks starting with an "l" are on the left side of the skull, 
landmarks starting with an "r" are on the right side. 

## Specifications for GDM analysis

When we do growth analyses, I like to put the older sample in the numerator spot and the younger sample in the denominator spot.  This means that the ratios reported are usually >1 which makes sense in terms of interpretation that things get larger with growth. 

So the way I set up the analysis in WinEDMA is:

- OLDER NUMERATOR FILE:     CZP0_mut_global.xyz
- YOUNGER NUMERATOR FILE:   CZEM_mut_global.xyz
- OLDER DENOMINATOR FILE:   CZP0_wt_global.xyz
- YOUNGER DENOMINATOR FILE: CZEM_wt_global.xyz


