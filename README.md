# GuamOxygenIsotopes
ISOLUTION and cave monitoring in Jinapsan Cave, Guam


MATLAB code in /ISOLUTION_MOD/ was originally written by Michael Deininger (2019):

It was then modified by Michael Deininger and Alexandra Noronha for this research.
NewDeininger, M., and Scholz, D. 2019. ISOLUTION 1.0: an ISOtope evoLUTION model describing the stable oxygen (δ18O) and carbon (δ13C) isotope values of speleothems. Int. J. of Speleology 48(1). 

New Files:
      BOUNDARY.m
            Additional model boundary values.
      RADIUS2.m
            Calculates calcite growth rate after Hansen et al., 2013, given cave T and PCO2, and Drip Ca.
      Run_Radius2_Script.m
            Runs RADIUS2 using the data given
      Run_Isotope_Calcite_Script.m
            Runs the ISOLUTION model from data files, bypassing user input steps.
Edited Files:
      CMODEL_FRAC.m
            Edited to add Affek and Zaarur 2014 calcite-water oxygen isotope fractionation.

References:

Affek, H.P., and Zaarur, S. 2014. Kinetic isotope effect in CO2 degassing: insight from clumped and oxygen isotopes in laboratory precipitation experiments. Geochimica et Cosmochimica Acta. 143, 319-330.	

Deininger, M., and Scholz, D. 2019. ISOLUTION 1.0: an ISOtope evoLUTION model describing the stable oxygen (δ18O) and carbon (δ13C) isotope values of speleothems. Int. J. of Speleology 48(1).

Hansen, M., Dreybrodt, W., and Scholz, D., 2013. Chemical evolution of dissolved inorganic carbon species flowing in thin water films and its implications for (rapid) degassing of CO2 during speleothem growth. Gheochim. Cosmochim. Acta 207, 242-251.
