function [I_ele, area_ele, areaS_ele] = sec_prop(ele_prop_ID,ele_num)

% The library for mechanical props of sections. If a new section needed,
% first add it here as a new case and then call it through the main file.

I_ele = zeros(length(ele_prop_ID),1);           % Moment of Inertia
area_ele = zeros(length(ele_prop_ID),1);        % Gross sectional area
areaS_ele = zeros(length(ele_prop_ID),1);       % Shear sectional area

j = 1;

for i = 1 : length(ele_prop_ID)
    switch ele_prop_ID(i)
        case 1                  % for W sections
            % w16x26
            I_ele(i) = 12528.566;
            area_ele(i) = 49.5483;
            areaS_ele(i) = 25.3064;
        case 2                  % for W sections
            % w16x31
            I_ele(i) = 15608.678;
            area_ele(i) = 58.8386;
            areaS_ele(i) = 28.1741;
        case 3                  % for circular sections with varying depth
            a = 4.7;
            b = 2.0;
            del = (a-b)/ele_num;
            step = del * j;
            rad = a - step + del/2;
            area_ele(i) = rad.^2*pi;            
            I_ele(i) = pi*rad.^4/4;   
            areaS_ele(i) = 0.9*area_ele(i);
            j = j + 1;
        case 4                  % for rectangular sections
            b = 5;
            h = 12;
            area_ele(i) = b*h;            % Cross sectional area
            I_ele(i) = b*h^3/12;          % Moment of Intertia
            areaS_ele(i) = 5*area_ele(i)/6;
        case 5                  % for circular section bars [I is trivial]
            area_ele(i) = 3;
            I_ele(i) = 1;
            areaS_ele(i) = 0.9*area_ele(i);
        case 6                  % for rectangular sections
            b = 20;
            h = 50;
            area_ele(i) = b*h;            % Cross sectional area
            I_ele(i) = b*h^3/12;          % Moment of Intertia
            areaS_ele(i) = 5*area_ele(i)/6;            
    end
end




