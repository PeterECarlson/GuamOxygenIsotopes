function output = ISOLUTION
disp('The ISOtope evoLUTION model - calculating the stable C and O isotope composition of calcite. Version 1.0 (18.02.2014, Michael Deininger)')
disp('Note that all inserts must be confirmed by pressing the ENTER botton.')
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
mode_to_analyse = input('If you want to analyse (1) a single isotope value, (2) the evolution of d18O and d13C in dependence on climatic parameter or (3) to simulate an artificial stalagmite record, please enter the related number and press ENTER.');
disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp(' ')
disp(' ')

if mode_to_analyse == 1
    mode = 'single';
elseif mode_to_analyse == 2
    mode = 'evolution';
elseif mode_to_analyse == 3
    mode = 'artificial';
else
    error('Wrong input value: try again!')
end

switch(mode)
    case 'single'
        disp('Please insert the following input parameter:')
        disp(' ')
        temperature = input('Temperature in °C: ');
        drip_interval = input('Drip interval in s: ');
        drip_ca = input('Drip water Ca in ppm: ');
        cave_pco2 = input('Cave air pCO2 in ppm: ');
        cave_pco2 = cave_pco2/1000000;
        h = input('Relative humidity (0<h<=1) (e.g. 90% = 0.9): ');
        while h >1 || h < 0
            h = input('h is not in the interval between 0 and 1, please insert the correct value for h: ');
        end
        if h == 1
            h = 0.999999;
        end
        v = input('Wind velocity inside the cave in m/s: ');
        phi = input('Mixing parameter phi (0<phi<=1) (for default value use 1): ');
        while phi >1 || phi < 0
            phi= input('Phi is not in the interval between 0 and 1, please insert the correct value for phi: ');
        end
        d18Oini = input('d18O value of the drip water (here of H2O) in per mil: ');
        d13Cini = input('d13C value of the drip water (here of HCO3-) in per mil: ');
        disp(' ')
        disp(' ')
        disp(' ')

      
        tmp_output = ISOTOPE_CALCITE(drip_interval, temperature, drip_ca,cave_pco2, h, v, phi, d18Oini, d13Cini);
        
        
        d18O_calcite = round(tmp_output(1)*100)/100;
        d13C_calcite = round(tmp_output(2)*100)/100;
        d18O_water = round(tmp_output(3)*100)/100;
        
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(sprintf('The d18O value of the precipitated calcite is %d per mil', d18O_calcite));
        disp(sprintf('The d13C value of the precipitated calcite is %d per mil', d13C_calcite));
        disp(sprintf('The d18O value of the water is %d per mil', d18O_water));
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(' ')
        disp(' ')
        disp(' ')
        
        output = [tmp_output(1), tmp_output(2), tmp_output(3)];
            
    case 'evolution'
        
        mode_for_evolution = input('If you want to analyse the evolution of d18O and d13C in dependence on (1) drip interval, (2) temperature, (3) pCO2cave, (4) pCO2dripwater or (5) the evaporation rate, ENTER the related number and press ENTER.');
        
        if mode_for_evolution == 1
            mode_evo = 'evo_dripinterval';
        elseif mode_for_evolution == 2
            mode_evo = 'evo_temperature';
        elseif mode_for_evolution == 3
            mode_evo = 'evo_pco2cave';
        elseif mode_for_evolution == 4
            mode_evo = 'evo_pco2dripwater';
        elseif mode_for_evolution == 5
            mode_evo = 'evo_evaporation';
        else
           error('Wrong input value: try again!');
        end
        
        switch(mode_evo)
            
            case 'evo_dripinterval'
                
                disp('Please insert the following input parameter:')
                disp(' ')
                disp('First, the drip interval range must be determined.')
                drip_start = input('Please insert the smallest value of the drip interval range; the smallest value possible is 1 second (the value must be an integer, i.e., no comma values like 10.5). E.g. 10 for 10 seconds');
                drip_end = input('Please insert the highest value of the drip interval range (e.g. 1800).');
                while drip_end < drip_start
                    input('The highest value of the drip interval should be higher then the smallest of the drip interval, please insert a higher value.')
                end
                drip_range = [drip_start:1:drip_end];
                disp(' ')
                disp(' ')

                disp('Please insert the following input parameter:')
                disp(' ')
                temperature = input('Temperature in °C, e.g., 10 for 10°C: ');
                drip_ca = input('Drip water Ca in ppm: ');
                cave_pco2 = input('Cave air pCO2 in ppm, e.g., 300 for 300ppm: ');
                cave_pco2 = cave_pco2/1000000;
                h = input('Relative humidity (0<h<=1) (e.g. 90% = 0.9): ');
                while h >1 || h < 0
                    h = input('h is not in the interval between 0 and 1, please insert the correct value for h: ');
                end
                if h == 1
                    h = 0.999999;
                end
                v = input('Wind velocity inside the cave in m/s: ');
                phi = input('Mixing parameter phi (0<phi<=1) (for default value use 1): ');
                while phi >1 || phi < 0
                    phi= input('Phi is not in the interval between 0 and 1, please insert the correct value for phi: ');
                end
                d18Oini = input('d18O value of the drip water (here of H2O) in per mil: ');
                d13Cini = input('d13C value of the drip water (here of HCO3-) in per mil: ');
                disp(' ')
                disp(' ')
                disp(' ')

                d18O_calcite = [];
                d13C_calcite = [];
                d18O_water = [];

                for j = 1:length(drip_range)
                    tmp_output = ISOTOPE_CALCITE(drip_range(j), temperature, drip_ca, cave_pco2, h, v, phi, d18Oini, d13Cini);
                    d18O_calcite(end+1) = tmp_output(1);
                    d13C_calcite(end+1) = tmp_output(2);
                    d18O_water(end+1) = tmp_output(3);
                end

                output = [d18O_calcite, d13C_calcite, d18O_water];

                figure(1), hold on, plot(drip_range, d18O_calcite,  'b', 'linewidth', 2), xlabel('drip interval (s)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(2), hold on, plot(drip_range, d13C_calcite, 'r', 'linewidth', 2), xlabel('drip interval (s)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{13}C_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(3), hold on, plot(drip_range, d18O_water, 'c', 'linewidth', 2), xlabel('drip interval (s)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{water} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)

            case 'evo_temperature'
                              
                disp('Please insert the following input parameter:')
                disp(' ')
                disp('First, the temperature range must be determined.')
                temperature_start = input('Please insert the smallest value of the temperature range (e.g. 5 for 5°C)');
                temperature_end = input('Please insert the highest value of the temperature range (e.g. 10 for 10°C)');
                while temperature_end < temperature_start
                    input('The highest value of the temperature range should be higher then the smallest temperature of the temperature range: please insert a higher value.')
                end
                temperature_range = [temperature_start:0.1:temperature_end];
                disp(' ')
                disp(' ')

                disp('Please insert the following input parameter:')
                disp(' ')
                dripinterval = input('Drip interval in seconds (no comma values), e.g., 300 for 300seconds: ');
                drip_pco2 = input('Drip water pCO2 in ppm, e.g., 10000 for 10.000ppm: ');
                drip_pco2 = drip_pco2/1000000;
                cave_pco2 = input('Cave air pCO2 in ppm, e.g., 300 for 300ppm: ');
                cave_pco2 = cave_pco2/1000000;
                h = input('Relative humidity (0<h<=1) (e.g. 90% = 0.9): ');
                while h >1 || h < 0
                    h = input('h is not in the interval between 0 and 1, please insert the correct value for h: ');
                end
                if h == 1
                    h = 0.999999;
                end
                v = input('Wind velocity inside the cave in m/s: ');
                phi = input('Mixing parameter phi (0<phi<=1) (for default value use 1): ');
                while phi >1 || phi < 0
                    phi= input('Phi is not in the interval between 0 and 1, please insert the correct value for phi: ');
                end
                d18Oini = input('d18O value of the drip water (here of H2O) in per mil: ');
                d13Cini = input('d13C value of the drip water (here of HCO3-) in per mil: ');
                disp(' ')
                disp(' ')
                disp(' ')

                d18O_calcite = [];
                d13C_calcite = [];
                d18O_water = [];

                for j = 1:length(temperature_range)
                    tmp_output = ISOTOPE_CALCITE(dripinterval, temperature_range(j), drip_pco2,cave_pco2, h, v, phi, d18Oini, d13Cini);
                    d18O_calcite(end+1) = tmp_output(1);
                    d13C_calcite(end+1) = tmp_output(2);
                    d18O_water(end+1) = tmp_output(3);
                end

                output = [d18O_calcite, d13C_calcite, d18O_water];

                figure(1), hold on, plot(temperature_range, d18O_calcite,  'b', 'linewidth', 2), xlabel('temperature (°C)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(2), hold on, plot(temperature_range, d13C_calcite, 'r', 'linewidth', 2), xlabel('temperature (°C)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{13}C_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(3), hold on, plot(temperature_range, d18O_water, 'c', 'linewidth', 2), xlabel('temperature (°C)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{water} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)

            case 'evo_pco2cave'
                
                disp('Please insert the following input parameter:')
                disp(' ')
                disp('First, the pCO2cave range must be determined.')
                cave_pco2_start = input('Please insert the smallest value of the cave air pCO2 (e.g. 300 for 300ppm)');
                cave_pco2_end = input('Please insert the highest value of cave air pCO2 (e.g. 5000 for 5000ppm)');
                while cave_pco2_end < cave_pco2_start
                    input('The highest value of the temperature range should be higher then the smallest temperature of the temperature range: please insert a higher value.')
                end
                pco2air_range = [cave_pco2_start:25:cave_pco2_end];
                pco2air_range = pco2air_range/1000000;
                disp(' ')
                disp(' ')
                disp('Please insert the following input parameter:')
                disp(' ')
                dripinterval = input('Drip interval in seconds (no comma values), e.g., 100 for 100seconds: ');
                temperature = input('Temperature in °C, e.g., 10 for 10°C: ');
                drip_pco2 = input('Drip water pCO2 in ppm, e.g., 10000 for 10.000ppm: ');
                drip_pco2 = drip_pco2/1000000;
                h = input('Relative humidity (0<h<=1) (e.g. 90% = 0.9): ');
                while h >1 || h < 0
                    h = input('h is not in the interval between 0 and 1, please insert the correct value for h: ');
                end
                if h == 1
                    h = 0.999999;
                end
                v = input('Wind velocity inside the cave in m/s: ');
                phi = input('Mixing parameter phi (0<phi<=1) (for default value use 1): ');
                while phi >1 || phi < 0
                    phi= input('Phi is not in the interval between 0 and 1, please insert the correct value for phi: ');
                end
                d18Oini = input('d18O value of the drip water (here of H2O) in per mil: ');
                d13Cini = input('d13C value of the drip water (here of HCO3-) in per mil: ');
                disp(' ')
                disp(' ')
                disp(' ')

                d18O_calcite = [];
                d13C_calcite = [];
                d18O_water = [];

                for j = 1:length(pco2air_range)
                    tmp_output = ISOTOPE_CALCITE(dripinterval, temperature, drip_pco2,pco2air_range(j), h, v, phi, d18Oini, d13Cini);
                    d18O_calcite(end+1) = tmp_output(1);
                    d13C_calcite(end+1) = tmp_output(2);
                    d18O_water(end+1) = tmp_output(3);
                end

                output = [d18O_calcite, d13C_calcite, d18O_water];

                figure(1), hold on, plot(pco2air_range*1000000, d18O_calcite,  'b', 'linewidth', 2), xlabel('cave air pCO_{2} (ppm)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(2), hold on, plot(pco2air_range*1000000, d13C_calcite, 'r', 'linewidth', 2), xlabel('cave air pCO_{2} (ppm)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{13}C_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(3), hold on, plot(pco2air_range*1000000, d18O_water, 'c', 'linewidth', 2), xlabel('cave air pCO_{2} (ppm)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{water} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)

            case 'evo_pco2dripwater'
            
                disp('Please insert the following input parameter:')
                disp(' ')
                disp('First, the pCO2dripwater range must be determined.')
                drip_pco2_start = input('Please insert the smallest value of the drip water pCO2 (e.g. 2000 for 2000ppm)');
                drip_pco2_end = input('Please insert the highest value of the drip water pCO2 (e.g. 10000 for 10.000ppm)');
                while drip_pco2_end < drip_pco2_start
                    input('The highest value of the temperature range should be higher then the smallest temperature of the temperature range: please insert a higher value.')
                end
                pco2drip_range = [drip_pco2_start:25:drip_pco2_end];
                pco2drip_range = pco2drip_range/1000000;
                disp(' ')
                disp(' ')
                disp('Please insert the following input parameter:')
                disp(' ')
                dripinterval = input('Drip interval in seconds (no comma values), e.g., 300 for 300seconds: ');
                temperature = input('Temperature in °C, e.g., 10 for 10°C: ');
                cave_pco2 = input('Cave air pCO2 in ppm, e.g., 300 for 300ppm: ');
                cave_pco2 = cave_pco2/1000000;
                h = input('Relative humidity (0<h<=1) (e.g. 90% = 0.9): ');
                while h >1 || h < 0
                    h = input('h is not in the interval between 0 and 1, please insert the correct value for h: ');
                end
                if h == 1
                    h = 0.999999;
                end
                v = input('Wind velocity inside the cave in m/s: ');
                phi = input('Mixing parameter phi (0<phi<=1) (for default value use 1): ');
                while phi >1 || phi < 0
                    phi= input('Phi is not in the interval between 0 and 1, please insert the correct value for phi: ');
                end
                d18Oini = input('d18O value of the drip water (here of H2O) in per mil: ');
                d13Cini = input('d13C value of the drip water (here of HCO3-) in per mil: ');
                disp(' ')
                disp(' ')
                disp(' ')

                d18O_calcite = [];
                d13C_calcite = [];
                d18O_water = [];

                for j = 1:length(pco2drip_range)
                    tmp_output = ISOTOPE_CALCITE(dripinterval, temperature, pco2drip_range(j), cave_pco2 ,h, v, phi, d18Oini, d13Cini);
                    d18O_calcite(end+1) = tmp_output(1);
                    d13C_calcite(end+1) = tmp_output(2);
                    d18O_water(end+1) = tmp_output(3);
                end

                output = [d18O_calcite, d13C_calcite, d18O_water];

                figure(1), hold on, plot(pco2drip_range*1000000, d18O_calcite,  'b', 'linewidth', 2), xlabel('drip water pCO_{2} (ppm)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(2), hold on, plot(pco2drip_range*1000000, d13C_calcite, 'r', 'linewidth', 2), xlabel('drip water pCO_{2} (ppm)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{13}C_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(3), hold on, plot(pco2drip_range*1000000, d18O_water, 'c', 'linewidth', 2), xlabel('drip water pCO_{2} (ppm)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{water} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
            
            case 'evo_evaporation'
                
                disp('Please insert the following input parameter:')
                disp(' ')
                disp('First, the relative humidity range must be determined.')
                h_start = input('Please insert the smallest value of the relative humidity (e.g. 0.5 for 50%)');
                h_end = input('Please insert the highest value of the drip water pCO2 (e.g. 1 for for 100%)');
                while h_end < h_start
                    input('The highest value of the relative humidity range should be higher then the smallest value for relative humidity: please insert a higher value.')
                end
                h_range = [h_start:0.05:h_end];
                disp(' ')
                disp(' ')
                disp('Please insert the following input parameter:')
                disp(' ')
                dripinterval = input('Drip interval in seconds (no comma values), e.g., 300 for 300seconds: ');
                temperature = input('Temperature in °C, e.g., 10 for 10°C: ');
                drip_pco2 = input('Drip water pCO2 in ppm, e.g., 10000 for 10.000ppm: ');
                drip_pco2 = drip_pco2/1000000;
                cave_pco2 = input('Cave air pCO2 in ppm, e.g., 300 for 300ppm: ');
                cave_pco2 = cave_pco2/1000000;
                v = input('Wind velocity inside the cave in m/s: ');
                phi = input('Mixing parameter phi (0<phi<=1) (for default value use 1): ');
                while phi >1 || phi < 0
                    phi= input('Phi is not in the interval between 0 and 1, please insert the correct value for phi: ');
                end
                d18Oini = input('d18O value of the drip water (here of H2O) in per mil: ');
                d13Cini = input('d13C value of the drip water (here of HCO3-) in per mil: ');
                disp(' ')
                disp(' ')
                disp(' ')

                d18O_calcite = [];
                d13C_calcite = [];
                d18O_water = [];

                for j = 1:length(h_range)
                    h = h_range(j);
                    if h == 1
                        h = 0.999999999999999;
                    end
                    tmp_output = ISOTOPE_CALCITE(dripinterval, temperature, drip_pco2, cave_pco2 ,h, v, phi, d18Oini, d13Cini);
                    d18O_calcite(end+1) = tmp_output(1);
                    d13C_calcite(end+1) = tmp_output(2);
                    d18O_water(end+1) = tmp_output(3);
                end

                output = [d18O_calcite, d13C_calcite, d18O_water];

                figure(1), hold on, plot(h_range*100, d18O_calcite,  'b', 'linewidth', 2), xlabel('relative humidity (%)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(2), hold on, plot(h_range*100, d13C_calcite, 'r', 'linewidth', 2), xlabel('relative humidity (%)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{13}C_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                figure(3), hold on, plot(h_range*100, d18O_water, 'c', 'linewidth', 2), xlabel('relative humidity (%)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{water} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
                
        end
        
    case 'artificial'
        
        raw_data = importdata('input.txt');
        parameter_evolution = raw_data.data;
    
        artificial_results = ones(length(parameter_evolution(:,1)),5);
        time_interval = [];
        d18O_depth = [];
        d13C_depth = [];
        radii = [];
        dfb = [];
        growthrate = [];
        
        for j = 1:length(parameter_evolution(:,1))
            
            if j == length(parameter_evolution(:,1))
                tmp_time_interval = [parameter_evolution(j,1):100:parameter_evolution(j,1)+1000];
                time_interval = [time_interval tmp_time_interval];
            else
                tmp_time_interval = [parameter_evolution(j,1):100:parameter_evolution(j+1,1)-100];
                time_interval = [time_interval tmp_time_interval];
            end
            
            dripinterval = parameter_evolution(j,2);
            temperature = parameter_evolution(j,3);
            drip_pco2 = parameter_evolution(j,4)/1000000;
            cave_pco2 = parameter_evolution(j,5)/1000000;
            d18Oini = parameter_evolution(j,6);
            d13Cini = parameter_evolution(j,7);
            
            tmp_output = ISOTOPE_CALCITE(dripinterval, temperature, drip_pco2, cave_pco2 ,0.999999999999, 0, 1, d18Oini, d13Cini);
            
            artificial_results(j,1) = tmp_output(1);
            artificial_results(j,2) = tmp_output(2);
            
            tmp_output2 = RADIUS2(dripinterval, temperature, 1, drip_pco2, cave_pco2);
            
            artificial_results(j,3) = tmp_output2(2);
            artificial_results(j,4) = tmp_output2(4);
            
            for k = 1:length(tmp_time_interval)
               if j == 1 && k == 1
                   dfb(end+1) = tmp_output2(2)*100;
               else
                   dfb(end+1) = dfb(end)+tmp_output2(2)*100;
               end
               radii(end+1) = tmp_output2(4);
               d18O_depth(end+1) = tmp_output(1);
               d13C_depth(end+1) = tmp_output(2);
               growthrate(end+1) = tmp_output2(2);
            end
            
        end
        
        age = parameter_evolution(end,1) - time_interval;
        
        figure(1), plot(age, dfb, 'k', 'linewidth', 2), xlabel('age (years)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('dfb (cm)', 'fontname', 'timesnewroman', 'fontsize', 16)
        figure(2), plot(age, growthrate*10, 'k', 'linewidth', 2), xlabel('age (years)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('growthrate (mm/a)', 'fontname', 'timesnewroman', 'fontsize', 16)
        figure(3), plot(radii,dfb, 'k', 'linewidth', 2), ylabel('dfb (cm)', 'fontname', 'timesnewroman', 'fontsize', 16), xlabel('radii (cm)', 'fontname', 'timesnewroman', 'fontsize', 16)
        figure(4), plot(d18O_depth,dfb, 'b', 'linewidth', 2), ylabel('dfb (cm)', 'fontname', 'timesnewroman', 'fontsize', 16), xlabel('\delta^{18}O_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
        figure(5), plot(d13C_depth, dfb, 'r', 'linewidth', 2), ylabel('dfb (cm)', 'fontname', 'timesnewroman', 'fontsize', 16), xlabel('\delta^{13}C_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
    	figure(6), plot(age, d18O_depth, 'b', 'linewidth', 2), ylabel('age (years)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{18}O_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
        figure(7), plot(age, d13C_depth, 'r', 'linewidth', 2), ylabel('age (years)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{13}C_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
        figure(8), plot(d18O_depth, d13C_depth, 'ko', 'linewidth', 2), xlabel('\delta^{18}O_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16), ylabel('\delta^{13}C_{calcite} (per mil)', 'fontname', 'timesnewroman', 'fontsize', 16)
        
end
