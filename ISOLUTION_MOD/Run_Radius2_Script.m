data = csvread('guam_growth_data.csv');
holder = [];
for i = 1:1:length(data)
    dripinterval = data(i,1);
    temperature = data(i,2);
    phi = 1;
    drip_ca = data(i,3);
    cave_pco2 = data(i,4);
    tmp = RADIUS2(dripinterval, temperature, phi, drip_ca, cave_pco2);
    holder = vertcat(holder, tmp);
end
results = holder;
