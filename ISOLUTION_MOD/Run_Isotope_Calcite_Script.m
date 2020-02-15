guam_data = csvread('guam_data.csv');
holder = [];
for i = 1:1:length(guam_data)
    tmp = ISOTOPE_CALCITE(guam_data(i,2), guam_data(i,3), guam_data(i,4), guam_data(i,5)/1e6, 0.9999999, 0.01, 1, guam_data(i,6), -12.5);
    holder = vertcat(holder, tmp);
end
results = horzcat(guam_data(:,1), holder(:,1));
csvwrite('\\Client\C$\Users\alexandramagana\Data Scripting\Python\2015 Guam Fractionation Paper\isolution_output_coplen.csv',results)
