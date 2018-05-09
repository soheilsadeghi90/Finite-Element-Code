function [I_ele, area_ele, areaS_ele] = sec_prop(ele_prop_ID,ele_num)

% The library for mechanical props of sections. If a new section needed,
% first add it here as a new case and then call it through the main file.

I_ele = zeros(length(ele_prop_ID),1);           % Moment of Inertia
area_ele = zeros(length(ele_prop_ID),1);        % Gross sectional area
areaS_ele = zeros(length(ele_prop_ID),1);       % Shear sectional area

j = 1;

for i = 1 : length(ele_prop_ID)
    switch ele_prop_ID(i)
        case 1                 
            I_ele(i) = 0.00213e8;
            area_ele(i) = 0.160e4;
            areaS_ele(i) = area_ele(i);
        case 2                  % for W sections
            I_ele(i) = 31183.117;
            area_ele(i) = 104.6;
            areaS_ele(i) = 32.0;
        case 3                  % for circular sections with varying depth
            I_ele(i) = 3001.6;
            area_ele(i) = 40.8;
            areaS_ele(i) = 12;
        case 4                  % for rectangular sections
            b = 11;
            h = 11;
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




