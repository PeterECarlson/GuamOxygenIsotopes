function output = ISOTOPE_CALCITE(d, TC, drip_ca, pCO2cave, h, V, phi, d18Oini, d13Cini);
%(08.12.2010/m)
%SOURCE CODE to calculate the isotopic compostion of oxygen 23.11.2010/m

%New Version, including the mass(volume)conservation of the first box, by a
%fixed volume of 0.1ml, by taking a part of the mass(volume) of the mobile
%box.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constants
IN = BOUNDARY;                                                                              %Boundary values
eva = EVAPORATION(TC, h, V);                                                                %Evaporationrate (mol/l)
FRACS = CMODEL_FRAC(TC);                                                                    %Fractionation factors
TK = 273.15 + TC;                                                                           %Absolute temperature (K)
Z = 0.0001/(1.188e-011 * TC^3 - 1.29e-011 * TC^2 + 7.875e-009 * TC + 4.844e-008);           %Tau precipitation (s); t=d/a (according to Baker 98)
alpha = (1.188e-011 * TC^3 - 1.29e-011 * TC^2 + 7.875e-009 * TC + 4.844e-008);
T = 125715.87302 - 16243.30688*TC + 1005.61111*TC.^2 - 32.71852*TC.^3 + 0.43333*TC.^4;      %Tau buffering, after Dreybrodt and Scholz (2010)

%Concentrations and mol mass
% outputsoil = KONSTANTEN(TC, pCO2);                                                          %Concentrations of the spezies in the solution, with respect to soil pCO2
outputcave = KONSTANTEN(TC, pCO2cave);                                                   %Concentrations of the spezies in the solution, with respect to cave pCO2        

% HCOSOIL = outputsoil{3}(3);                                                                 %HCO3- concentration, with respect to soil pCO2 (mol/l)    
HCOCAVE = outputcave{3}(3)/sqrt(0.8);                                                                 %HCO3- concentration, with respect to cave pCO2 (mol/l)

h2o_ini = 0.1/18;                                                                           %Mol mass of the water, with respect to the volume of a single box (mol)
hco_ini = (drip_ca) * 2 * 1e-4;                                                                     %Mol mass of the HCO3-(soil), with respect to the volume of a single box (mol)
hco_eq = HCOCAVE*1e-4;                                                                      %Mol mass of the HCO3-(cave), with respect to the volume of a single box (mol)

%Fractionation facors for oxygen isotope
eps_m = FRACS.a18_m - 1;                                                                    %
avl = (-7356/TK + 15.38)/1000 + 1;
abl = 1/(FRACS.e18_hco_h2o + 1);
a = 1/1.008;
f = 1/6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rdrop18_w = ( (d18Oini / 1000) + 1) * IN.R18smow;
Rdrop18_b = Rdrop18_w / (FRACS.e18_hco_h2o + 1);
Rv18 = avl * Rdrop18_w;
Rdrop13_b = ( (d13Cini / 1000) + 1) * IN.R13vpdb;

% Definition of the input-parameter of the program "O18EVA" at the time
% t=0s, e.g. at the beginning of the mix process

hco_mix = hco_ini;
h2o_mix = h2o_ini;
HCOMIX = hco_mix / 1e-4;

r_hco18_mix = Rdrop18_b;
r_hco13_mix = Rdrop13_b;
r_h2o18_mix = Rdrop18_w;

number = 0;

r18res = 0;
r18mix = 1;

r13res = 0;
r13mix = 1;

% while number ~= 75
while r18mix ~= r18res && r13mix ~= r13res;
    
    number = number + 1;
    
    r18_hco_res = r_hco18_mix;
    r13_hco_res = r_hco13_mix;
    
    temp = O18EVA(d, TC, pCO2cave,  h, V, r_hco18_mix, r_h2o18_mix, Rv18, HCOMIX, h2o_mix, r_hco13_mix);

    hco_out = temp(4);      %mol mass of hco
    h2o_out = temp(6);      %mol mass of h2o
    
    r_hco18_out = temp(1);        %oxygen isotopic ratio of hco
    r_hco13_out = temp(7);        %carbon isotopic ratio of hco
    r_h2o18_out = temp(2);        %oxygen isotopic ratio of h2o
    
    %%% 1) simple mixprocess
    hco_mix = phi*hco_ini + (1-phi)*hco_out;        %mixing of hco
    h2o_mix = phi*h2o_ini + (1-phi)*h2o_out;        %mixing of h2o
    
    phi_r_b = 1 / (1 + (1-phi) / phi*hco_out/hco_ini);      %isotopic mixing parameter of hco
    phi_r_w = 1 / (1 + (1-phi) / phi*h2o_out/h2o_ini);      %isotopic mixing parameter of h2o
    
    r_hco18_mix = phi_r_b*Rdrop18_b + (1-phi_r_b)*r_hco18_out;        %new oxygen isotopic ratio of hco
    r_h2o18_mix = phi_r_w*Rdrop18_w + (1-phi_r_w)*r_h2o18_out;        %new oxygen isotopic ratio of h2o
    r_hco13_mix = phi_r_b*Rdrop13_b + (1-phi_r_b)*r_hco13_out;        %new carbon isotopic ratio of hco
    
    H2O_mix = h2o_mix*18/1000;  %mol -> l (mol*(g/mol)/(g/kg)*(l/kg)
    HCOMIX = hco_mix / H2O_mix;     %new hco condentration
    
    %end of the extended mixprocess
    
    r18mix = round(r_hco18_mix*10^13);
    r18res = round(r18_hco_res*10^13);
    
    r13mix = round(r_hco13_mix*10^13);
    r13res = round(r13_hco_res*10^13);
end

temp = O18EVA_MEAN(d, TC, pCO2cave, h, V, r_hco18_mix, r_h2o18_mix, Rv18, HCOMIX, h2o_mix, r_hco13_mix);

r_hco18 = temp(:,1);
r_hco13 = temp(:,5);
hco = temp(:,3);

delta_1 = temp(:,6);
delta_0 = delta_1(1);
delta_end = delta_1(end);

r_h2o18 = temp(:,2);
h2o = temp(:,4);

% TK = TC+273.16;
% e18_h2o_caco3 = (18030./TK - 32.42)/1000;
tmp = isnan(r_hco18);
if tmp ~= 1

r_hco18_mean = sum( (r_hco18(1:end-1) + (diff(r_hco18)/2)) .* ( (-diff(hco)) ./ (hco(1)-hco(end)) ) );
% d18Ocalcite = (((r_hco18_mean*(FRACS.e18_hco_h2o+1))*(e18_h2o_caco3+1))/IN.R18vpdb - 1)*1000;
d18Ocalcite = (r_hco18_mean*(FRACS.e18_hco_caco + 1)/IN.R18vpdb - 1)*1000;

r_hco13_mean = sum( (r_hco13(1:end-1) + (diff(r_hco13)/2)) .* ( (-diff(hco)) ./ (hco(1)-hco(end)) ) );
d13Ccalcite = (r_hco13_mean*(FRACS.e13_hco_caco + 1)/IN.R13vpdb - 1)*1000;

r_h2o18_mean = sum( (r_h2o18(1:end-1) + (diff(r_h2o18)/2)) .* ( (-diff(h2o)) ./ (h2o(1)-h2o(end)) ) );
d18Owater = (r_h2o18_mean/IN.R18smow - 1)*1000;

d18Ovapor = (Rv18/IN.R18smow - 1)*1000;

output = [d18Ocalcite, d13Ccalcite, d18Owater, d18Ovapor, hco_mix, h2o_mix, r_hco18_mean, r_hco13_mean, r_h2o18_mean, delta_0, delta_end];

else
output = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN];    
end