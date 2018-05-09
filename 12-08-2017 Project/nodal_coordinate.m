function [tag, x, y, p, u0, u_dot0] = nodal_coordinate()

% This function is one of the INPUT functions. Operator is responsible to fill it.
% Enter x and y and p. For five verifying examples, the inputs are
% prepared. Just uncomment if needed. 

% UNITS: cm and N are preferred. 

% --------------------- FIRST EXAMPLE -----------------------

x = [0  75  150 ];
y = [0  0   0   ];

p = [0    0   0
     0    0   0
     0    0   0];      % Px, Py and moment at node 1, respectively, in N and cm

% if dynamic analysis needed, enter the following initial conditions
u0 = [0 0   0;-0.0000 1.9632   0.0455;-0.0000 6.0000   0.0580];  % First Modal Shape
u_dot0 = [0 0   0;0 0   0;0 0   0];

% ------------------ SECOND EXAMPLE -------------------------
% ele_num = 5/0.1;
% del = 10;
% i = 0:del:300;
% y = [i 300*ones(1,20)];
% j = 0:del:200;
% x = [zeros(1,30) j(1:end)];
% p = [zeros(ele_num,3);100e3 0 0];      % Px, Py and moment at node 1, respectively, in N and cm
% u0 = zeros(ele_num+1,3);               
% u_dot0 = zeros(ele_num+1,3);

% ---------------------THIRD EXAMPLE -----------------------

% x = [0 0 0 0 0 0 0 600 600 600 600 600 600 600 200 400 200 400 200 400];
% y = [0 200 400 560 720 880 1040 0 200 400 560 720 880 1040 400 400 720 720 1040 1040];
% p = zeros(20,3); % Px, Py and moment at node 1, respectively, in N and cm
% p(7,1) = 7000;
% u0 = zeros(20,3);
% u_dot0 = zeros(20,3);

% --------------------- FOURTH EXAMPLE -------------------------
% x = [0 50 100 150 200 250 300 350 400 200 200 200 200];
% y = [0 0 0 0 0 0 0 0 0 -50 -100 -150 -200];
% 
% p = [0      10000   0
%      0      0       0
%      0      0       0
%      0      0       0
%      0      0       0
%      0      0       0
%      0      0       0
%      0      0       0
%      0      0       0
%      0      0       0
%      0      0       0
%      0      0       0
%      0      0       0];      % Px, Py and moment at node 1, respectively, in N and cm
%  
% u0 = zeros(13,3);
% u_dot0 = zeros(13,3);

 % --------------------- FIFTH EXAMPLE -------------------------

% x = [0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500];
% y = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
% 
% p = zeros(26,3);      % Px, Py and moment at node 1, respectively, in N and cm
% p([8,19],:) = [0 20000 0; 0 -20000 0];
% u0 = zeros(26,3);
% u_dot0 = zeros(26,3);

% -------------------------------------------------------------
 n = length(x);
 tag = 1 : n;
 
end

