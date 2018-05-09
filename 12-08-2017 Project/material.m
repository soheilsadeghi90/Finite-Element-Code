function [MAT] = material(MAT_ID)

% Material library. Each row in the MAT matrix consists three values. First
% is Elastic modulus, the second one is nu and the last one is the density. MAT = [E, v, ro]

% Initialization
MAT = zeros(length(MAT_ID),3);

% insert values in cm and N
for i = 1:length(MAT_ID)
    switch MAT_ID(i)
        case 1
            MAT(i,1) = 20000000;       % Elastic Modulus (N/cm2)
            MAT(i,2) = 0.29;
            MAT(i,3) = 7.80e-5;       % gram/cm3 divided by 100 to convert to N/cm2
       case 2
            MAT(i,1) = 20000000;
            MAT(i,2) = 0.29;
            MAT(i,3) = 7.86e-5;       % gram/cm3 divided by 100 to convert to N/cm2
       case 3
            MAT(i,1) = 21000000;
            MAT(i,2) = 0.365;     
            MAT(i,3) = 7.85e-5;       % gram/cm3 divided by 100 to convert to N/cm2
       case 4
            MAT(i,1) = 20000000;
            MAT(i,2) = 0.3;
            MAT(i,3) = 7.85e-5;       % gram/cm3 divided by 100 to convert to N/cm2            
    end
end

end
