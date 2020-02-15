function output = BOUNDARY
%(08.12.2010/m)
%constant parameter

pCO2soil = 10000e-6;
pCO2cave = 600e-6;

R2smow = 0.00015575;
R13vpdb = 0.0112372;
R18smow = 0.0020052;
R18vpdb = 0.0020672;

output = struct('pCO2soil', pCO2soil, 'pCO2cave', pCO2cave, 'R13vpdb', R13vpdb, 'R18smow', R18smow, 'R18vpdb', R18vpdb, 'R2smow', R2smow);
