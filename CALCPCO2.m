function output = CALCPCO2(ca, TC)

%input variable is Ca2+ concentration of the drip water in mol/l and the
%ave air temperature in °C

pco2int = [0:100000e-6:1000000e-6];

run = 1;    %varialbe to stop the while loop when set to 0

while run==1
    
    ca_int = [];
    
    for i = 1:length(pco2int)
    
         concentrations = KONSTANTEN(TC, pco2int(i));
         ca_int(end+1) = concentrations{3}(1);
    
    end
    
    [next_pco2, index_pco2] = min(abs(ca_int-ca));
    
    if ca_int(index_pco2)<ca
        pco2int = [pco2int(index_pco2):(pco2int(index_pco2+1)-pco2int(index_pco2))/10:pco2int(index_pco2+1)];
    else
        pco2int = [pco2int(index_pco2-1):(pco2int(index_pco2)-pco2int(index_pco2-1))/10:pco2int(index_pco2)]; 
    end
    
    if ca_int(index_pco2)==ca
        run=0;
    end
    
end

output = pco2int(index_pco2)*1e6;

end