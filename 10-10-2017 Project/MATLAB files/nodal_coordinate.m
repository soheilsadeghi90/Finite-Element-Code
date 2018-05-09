function [tag, x, y, p] = nodal_coordinate()

% This function is one of the INPUT functions. Operator is responsible to fill it.
% Enter x and y and p. For five verifying examples, the inputs are
% prepared. Just uncomment if needed. 

% UNITS: cm and N are preferred. 

% ------------------ FIRST EXAMPLE -------------------------
% ele_num = 50;
% j = 1:ele_num;
% del = 350/ele_num;
% j = del*j;
% x = [0  j];
% y = zeros(1,ele_num+1);
% p = [zeros(ele_num,3);6*10^5 0 0];      % Px, Py and moment at node 1, respectively, in N and cm

% ---------------------SECOND EXAMPLE -----------------------

% L = 300;
% x = [0  0   L];
% y = [0  280 280];
% p = [0    0   0;0 0   0;0 -2000   0];      % Px, Py and moment at node 1, respectively, in N and cm

% --------------------- THIRD EXAMPLE -------------------------

% x = [0	700 0   700 0   700 0   700];
% y = [0	0   300 300 600 600 900 900];
% 
% p = [0      0   0
%      0      0   0
%      10000  0   0
%      0      0   0
%      20000  0   0
%      0      0   0
%      30000  0   0
%      0      0   0];      % Px, Py and moment at node 1, respectively, in N and cm

 % --------------------- FOURTH EXAMPLE -------------------------

% x = [0  400 800 1200    1600    1200    800 400];
% y = [0  0   0   0       0       300     300 300];
% 
% p =[    0   0       0   
%         0   -3000   0
%         0   -4000   0
%         0   -3000   0
%         0   0       0
%         0   0       0
%         0   0       0
%         0   0       0       ];      % Px, Py and moment at node 1, respectively, in N and cm

% --------------------- FIFTH EXAMPLE -------------------------

x = [0  500 750 1000    1500 ];
y = [0  0   0   0       0    ];

p = [    0   0        0
        0   0        50000
        0   -60000   0
        0   0        0
        0   0        0    ];          % Px, Py and moment at node 1, respectively, in N and cm


% --------------------- SIXTH EXAMPLE -------------------------

% x = [0  600 1200];
% y = [0  0   0];
% 
% p = [    0   0        0
%          0   0        0
%          0   0   0 ];          % Px, Py and moment at node 1, respectively, in N and cm

% -------------------------------------------------------------
 n = length(x);
 tag = 1 : n;
 
end

