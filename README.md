# AAT_SEP_sourceselection
These codes can be used to select science targets, sky fibre positions, calibration and guide stars for the AAT field selection. 

Varsrcprep generates the science targets.

PutsrcinSEPreg is for adding a particular set of other sources to some previously selected fields. 

SelectWDs1 is for selecting white dwarfs for the spectral calibration. 

Guidestargen is for generating a list of guide stars. For this, you first have to generate a list of sources from SIMBAD within 1 degree of the selected field center, and get it to list the source name, RA, Dec, r mag, and proper motions. 

skyfibrefromJacob generates the sky fibre positions from the fits file that Jacob Ider Chitham created, which is also attached. 

testplotfldfiles runs some checks on the final fld files, to check for duplicates, and if any selected source lies outside of the selected fields. 
