function output = O18EVA(tmax, TC, pCO2cave, h, V, R18_hco_ini, R18_h2o_ini, R18v, HCOMIX, h2o_new, R13_hco_ini);

% Sourcecode to develope the evolution of the isotopicratio of the oxygen
% compostion of the oxygen isotopes 16O and 18O as a function of
% temperature TC, supersaturation (pCO2), relative humidity (h) and wind
% velocity (v). %(08.12.2010/m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Constants
IN = BOUNDARY;                                                                              %Boundary values
eva = EVAPORATION(TC, h, V);                                                                %Evaporationrate (mol/l)
FRACS = CMODEL_FRAC(TC);                                                                    %Fractionation factors
TK = 273.15 + TC;                                                                           %Absolute temperature (K)
alpha_p = (1.188e-011 * TC^3 - 1.29e-011 * TC^2 + 7.875e-009 * TC + 4.844e-008);           %Tau precipitation (s); t=d/a (according to Baker 98)
T = 125715.87302 - 16243.30688*TC + 1005.61111*TC.^2 - 32.71852*TC.^3 + 0.43333*TC.^4;      %Tau buffering, after Dreybrodt and Scholz (2010)

%Concentrations and mol mass
outputcave = KONSTANTEN(TC, pCO2cave);                                                   %Concentrations of the spezies in the solution, with respect to cave pCO2                            
HCOCAVE = outputcave{3}(3)/sqrt(0.8);                                                                 %HCO3- concentration, with respect to cave pCO2 (mol/l)

h2o_ini = h2o_new;                                                                           %Mol mass of the water, with respect to the volume of a single box (mol)
H20_ini = h2o_ini*18/1000;
hco_ini = HCOMIX*H20_ini;                                                                      %Mol mass of the HCO3-(soil), with respect to the volume of a single box (mol)
hco_eq = HCOCAVE*H20_ini;                                                                      %Mol mass of the HCO3-(cave), with respect to the volume of a single box (mol)

%Fractionation facors for oxygen isotope
eps_m = FRACS.a18_m - 1;                                                                    %
avl = ((-7356/TK + 15.38)/1000 + 1);
abl = 1/(FRACS.e18_hco_h2o + 1);
a = 1/1.008*1.003;
f = 1/6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculation of the 18R

r_hco18 = R18_hco_ini;
r_h2o18 = R18_h2o_ini;
r_hco13 = R13_hco_ini;

if tmax > floor(h2o_ini/eva)
    tmax = floor(h2o_ini/eva);
    disp('DRIPINTERVALL IS TO LONG, THE WATERLAYER EVAPORATES COMPLETLY FOR THE GIVEN d')
end
    
dt = 1;
t=1:tmax;
        
HCO = [];                               %Konzentration von HCO3-
HCO(end+1) = HCOMIX;                    %Start-Konzentration der Lösung ohne Verdunstung
hco = [];                               %Menge an HCO3-
hco(end+1) = hco_ini;
H2O = [];
H2O(end+1) = H20_ini;
h2o = [];
h2o(end+1) = h2o_new;

%Die Gleichgewichtskonzentration ist unabhängig von der
%"Restwassermenge" und daher konstant
HCO_EQ = HCOCAVE;                                           %Equilibriumconcentration

for i=1:length(t)
    
    delta = (H2O(end)/1000)/0.001;
    
    %Neue Wassermenge
    h2o(end+1) = h2o(end) - eva*dt;                                   %Water (mol)
    d_h2o = -eva;                                               %Evaporationrate (mol/l)
    H2O(end+1) = h2o(end)*18*1e-3;                                   %Water (l)
            
    %Berechnung der Konzentration unter Berücksichtigung der
    %Verdundstung
    HCO_temp = (HCO(end) - HCO_EQ) * exp(-dt/(delta/alpha_p)) + HCO_EQ;       %HCO3- concentration after timeintervall dt
    hco(end+1) = HCO_temp * H2O(end-1);                         %HCO3- mass (mol)
            
    HCO(end+1) = HCO_temp * (H2O(end-1)/H2O(end));              %HCO3- concentration after timeintervall dt and the evaporation of water
            
    r_hco18(end+1,1) = r_hco18(end) + ((eps_m*(hco(end)-hco(end-1))/hco(end)-1/T) * r_hco18(end) + abl/T*r_h2o18(end)) * dt;
    r_h2o18(end+1,1) = r_h2o18(end) + ((hco(end)/h2o(end)/T - f/abl/h2o(end)*(hco(end)-hco(end-1))) * r_hco18(end) + (d_h2o/h2o(end)*(a*avl/(1-h)-1) - hco(end)/h2o(end)*abl/T) * r_h2o18(end) - a*h/(1-h)*R18v/h2o(end)*d_h2o) * dt;
    r_hco13(end+1,1) = r_hco13(end) + ((FRACS.a13_m - 1)*(hco(end)-hco(end-1))/hco(end)) * r_hco13(end);
                        
end

output = [r_hco18(end), r_h2o18(end), HCO(end), hco(end), H2O(end), h2o(end), r_hco13(end)];
