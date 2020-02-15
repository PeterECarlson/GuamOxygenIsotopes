function output = CMODEL_FRAC(TC)
%fractionation and enrichmant factors (08.12.2010/m)

TK = TC + 273.16;                                                            %absolute temperature (K)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Version using epsilon = exp(...)-1 = alpha - 1 = epsilon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%carbon C (absolute values)
e13_b3 = exp((- 9483/TK + 23.89)/1000)-1;                                   %Mook 1974 %HCO3 -> CO2^{g}            (1000ln alpha)
e13_d3 = ((- 4232./TK + 15.10) / 1000);                                     %Mook 2000 %HCO3 -> CaCO3
e13_m = 0.5 * (e13_b3 + e13_d3);                                            %mean enrichment factor for Rayleighfractionation HCO3- -> CO2 + CaCO3
a13_m = e13_m + 1;                                                          %mean enrichment factor for Rayleighfractionation HCO3- -> CO2 + CaCO3

%oxygen O (absolute values)
e18_c4 = exp((- 2590000./(TK.^2) - 1.89)/1000)-1;                           %Beck 2005 %HCO3 -> H2O                 (1000ln alpha)
e18_b5 = exp((2520000./(TK.^2) + 12.12)/1000)-1;                            %Beck 2005 %H2O -> CO_2^{aq}            (1000ln alpha) 
e18_a1 = (-160515./(TK.^2)+1441.76./TK - 1.9585)/1000;                      %Thorstenson 2004 %CO_2^{g} -> CO_2^{aq}


e18_s1 = ((e18_c4+1)*(e18_b5+1)/(e18_a1+1))-1;                              %[Thorstenson 2004; Beck 2005] %HCO3 -> CO_2{g}  

% e18_d1 = (2780000./(TK.^2) - 3.39) / 1000;                                %O'Neil 1969 %H2O -> CaCO3
% e18_d1 = exp(((15630./TK) - 23.29) /1000)-1;                              %Affek  and Zaarur 2014 %H20 -> CaCO3   (1000ln alpha)
e18_d1 = exp((17400./TK - 28.6)/1000)-1;                                    %Coplen 2007
% e18_d1 = exp((18030./TK - 32.42)/1000)-1;                                 %Kim and O'Neil 1997 %H2O -> CaCO3      (1000ln alpha)


e18_s2 = ((e18_c4+1)*(e18_d1+1))-1;                                         %[Coplen 2007; Beck 2005] %HCO3 -> CaCO3

e18_m = 2/6 * e18_s1 + 3/6 * e18_s2 + 1/6 * e18_c4;                         %[Thorstenson 2004; Beck 2005; Coplen 2007] %HCO3- -> CO2 + CaCO3 + H2O %mean enrichment for Rayleighfractionation
a18_m = e18_m + 1;                                                          %[Thorstenson 2004; Beck 2005; Coplen 2007] %HCO3- -> CO2 + CaCO3 + H2O %mean fractionationfactor for Rayleighfractionation

%Output: fractionation and enrichment factors
output = struct('e13_hco_caco', e13_d3, 'a13_m', a13_m, 'e18_hco_caco', e18_s2, 'e18_hco_h2o', e18_c4, 'a18_m', a18_m);
