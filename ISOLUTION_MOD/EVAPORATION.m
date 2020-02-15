function output = EVAPORATION(TC, h, V)
%calculation of the evaporationrate, Bansall and Xie 1998, (08.12.2010/m)

T = 273.15 + TC;                                                                                %absolute temperature (K)
E = exp(-6094.4642/T + 21.1249952 - 2.724552e-2*T + 1.6853396e-5*T^2 + 2.4575506*log(T));       %D. Sonntag (1982) [Pa]
E = E/1.3332e2;                                                                                 %Pascal->Torr Stöcker 2005 [mmHg ]
A = pi*0.017841241^2;                                                                           %Oberfläche einer Box mit Volumen von 0.1ml und einer Schichtdicke von 0.1mm

%eva = A * (0.002198 + 0.0398*0^0.5756) * E/1.3332e2 * (1-h)*1e6/3600; %µl/s (Bansal1998)
eva = A * (0.002198 + 0.0398*V^0.5756) * E * (1-h)*1000/3600/18; %mol/s (Bansal1998)

%Output: evaporationrate (mol/s)
output = eva;
