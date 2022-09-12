

2022-08-07 AB This repository includes the codes used in Bonsor et al, 2022, Nature Astronomy (please note the relationship between the collision code used here and in Bonsor et al, 2023, MNRAS, in prep)

Fig 1: 
bonsor_2022_cf_twd.ipynb
Plots two samples of white dwarfs (known as Sample One and Sample Two in the text, noting that there is overlap between them), showing their Ca/Fe as a function of the WD temperature. 
Bonsor_source_data_table1_all.csv: the combined sample teff, log (Ca/Fe), 1-sigma error in log (Ca/Fe)
bonsor_source_data_table1_core_mantle_rich.csv
Highlighted are those white dwarfs where all the elemental abundances considered together indicate the accretion of core or mantle--rich material as the most likely explanation, to >3\sigma over a standard model. 





Fig 4 
collevol_origin.c

This code models the collisional evolution of a planetesimal belt, containing particles from size D_min up to D_max, split into logarithmically spaced mass bins. It follows the rate of catastrophic and/or resurfacing collisions at each timestep, redistributing the mass accordingly.

collevol_origin.c : tracks where the material 'originated', as well as the collisional evolution. This allows you to answer the question (and similar): at a given timestep how much mass in e.g. 1km bodies are collision fragments of 100km bodies. [This version cannot also track volatiles.]

In order to run the code, you need to fill in the input parameters in inparam.in. 

gcc -lm -o collevol_origin collevol_origin.c
./collevol_origin >out & 

Output Files: 

collouts.dat 
colloutv.dat (for volatiles only)
origin****.dat (**** labels each diameter bin)

Fig 4 a plots the fraction of 30km planetesimals that are fragments of Plutos (D>1400km), for a belt at 1au with initial total mass in 100um to 5,000km particles of 100Mearth, with <e>=0.1. 



