

2022-08-07 AB This repository includes the codes used in Bonsor et al, 2022, Nature Astronomy and Bonsor et al, 2023, MNRAS, in prep

This code models the collisional evolution of a planetesimal belt, containing particles from size D_min up to D_max, split into logarithmically spaced mass bins. It follows the rate of catastrophic and/or resurfacing collisions at each timestep, redistributing the mass accordingly. There are two main versions:

collevol_origin.c : tracks where the material 'originated', as well as the collisional evolution. This allows you to answer the question (and similar): at a given timestep how much mass in e.g. 1km bodies are collision fragments of 100km bodies. [This version cannot also track volatiles.]

collevol_resurf_chid.c: tracks the collisional evolution of bodies made from both volatiles and solids. Volatiles are released to gas following each catastrophic and/or resurfacing collision. 

In order to run the code, you need to fill in the input parameters in inparam.in. [Note that at present the format is slightly different for collevol_origin.c and collevol_resurf_chi_d.c] 

gcc -lm -o collevol_origin collevol_origin.c
./collevol_origin >out & 

Output Files: 

collouts.dat 
colloutv.dat (for volatiles only)

