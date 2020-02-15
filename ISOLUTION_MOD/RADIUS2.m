function output = RADIUS2(dripinterval, temperature, phi, drip_ca, cave_pco2)
%Berechnung des Gleichgewichtsradius und der Wachstumsrate f�r einen Ca
%Excess von CASOIL-CACAVE unter Ber�cksichtigung von Mischungsprozessen
%(Anhang "Mix")
%(08.12.2010/m)

dT = dripinterval;
TC = temperature;

tmp_cave = KONSTANTEN(TC, cave_pco2*1e-6);

CASOIL = (drip_ca/(1e3 * 40.08));
CACAVE = tmp_cave{3}(1);

%(mol/L - mol/L) * 1 L/0.001 m^3 = mol/m^3
CaEx = ( CASOIL - CACAVE ) * 1e3;       %Kalziumexcess [mol/m�]

%m, film thickness
delta = 1e-5;

%m, drop volume
drop_volume = 5e-8;

%m/s, Equation from Hansen2013 pg 244
alpha = (0.52 + 0.04 * TC + 0.004 * TC^2) * 1e-7;

% m/(m/s) = s
Z = delta / alpha;

%g/mol * m/s * (1 - exp(-s/s)) * mol/m^3 = g/(m^2*s), Growth rate
W0 = 100.09 * delta / dT * ( 1 - exp(-dT / Z) ) * CaEx;

%h/day * min/day * s/day = s
seconds_per_day = 24*60*60;

%m^3/m = m^2
area = drop_volume/delta;

%g/(m^2*s) * m^2 * s/day = g/day
output = [W0*area*seconds_per_day];
