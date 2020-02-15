# GuamOxygenIsotopes
ISOLUTION and cave monitoring in Jinapsan Cave, Guam
      oxygen_isotope_stats_functions.py:

       
      oxygen_isotope_stats_functions.py:
            This contains functionationality to:
                  read the supplementary material
                  propagate uncertainties for some derived variables
                  correlate data using a 2-part piecewise linear function
                  calculate linear correlation significance, correction for autoregressive
                  characteristics (After Hu et al., 2017)
      
      References:
      
      Hu, J., Emile-Geay, J., Partin, J. 2017. Correlation-based interpretations of paleoclimate data—
            where statistics meet past climates. Earth and Plan. Sci. Let. 459:362-371.

      /ISOLUTION_MOD/*
            MATLAB code in /ISOLUTION_MOD/ was originally written by Michael Deininger (Deininger and Scholz, 2019). 
            It was then modified by Michael Deininger and Alexandra Noronha for this research. 

      New ISOLUTION Files:
            BOUNDARY.m
                  Additional model boundary values.
            RADIUS2.m
                  Calculates calcite growth rate after Hansen et al., 2013, given cave T and PCO2, and Drip Ca.
            Run_Radius2_Script.m
                  Runs RADIUS2 using the data given
            Run_Isotope_Calcite_Script.m
                  Runs the ISOLUTION model from data files, bypassing user input steps.
      Edited ISOLUTION Files:
            CMODEL_FRAC.m
                  Edited to add Affek and Zaarur 2014 calcite-water oxygen isotope fractionation.

      References:

      Affek, H.P., and Zaarur, S. 2014. Kinetic isotope effect in CO2 degassing: insight from clumped and
            oxygen isotopes in laboratory precipitation experiments. Geochimica et Cosmochimica Acta. 143, 319-330.	

      Deininger, M., and Scholz, D. 2019. ISOLUTION 1.0: an ISOtope evoLUTION model describing the stable oxygen 
            (δ18O) and carbon (δ13C) isotope values of speleothems. Int. J. of Speleology 48(1).

      Hansen, M., Dreybrodt, W., and Scholz, D., 2013. Chemical evolution of dissolved inorganic carbon species 
            flowing in thin water films and its implications for (rapid) degassing of CO2 during speleothem growth. 
            Geochim. Cosmochim. Acta 207, 242-251.
      
   
