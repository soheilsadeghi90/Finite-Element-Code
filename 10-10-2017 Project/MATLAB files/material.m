function [MAT] = material(MAT_ID)

% Material library. Each row in the MAT matrix consists two values. First
% is Elastic modulus and the second one is nu. MAT = [E, v]

% Initialization
MAT = zeros(length(MAT_ID),2);

% insert values in cm and N
for i = 1:length(MAT_ID)
    switch MAT_ID(i)
        case 1
            MAT(i,1) = 20000000;       % Elastic Modulus (N/cm2)
            MAT(i,2) = 0.3;
       case 2
            MAT(i,1) = 6800000;
            MAT(i,2) =0.3;
       case 3
            MAT(i,1) = 21000000;
            MAT(i,2) =0.29;     
       case 4
            MAT(i,1) = 29000;
            MAT(i,2) =0.3;                 
    end
end

end
