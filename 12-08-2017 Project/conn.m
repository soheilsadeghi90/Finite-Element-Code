function [conn,supp,nodal_disp] = conn()

% This file is an INPUT file. Operator is responsible to enter inputs in
% this file. For five verifying examples, inputs are prepared. Uncomment
% each if needed.

% conn: connectivity matrix. Each row includes element tag, xconnection,
% yconnection, rotational connection for the starting node and same
% connections for the ending nodes. So each rows has to have 7 inputs. [if
% restrained 1, if not 0]

% supp: nodal support. Each row includes nodal tag, x constraint, y
% constraint and rotational constraints. So each row has to have 4 inputs. [if
% restrained 1, if not 0]

% nodal_disp: in nodal displacement in applied at each point, it can be
% entered through this matrix. Each row has a node tag, x disp, y disp and
% rotation. So each row has 4 inputs. 

% ------------------- FIRST EXAMPLE ---------------------- 

conn = [1   1   1   1   1   1   1
        1   1   1   1   1   1   1];
    
supp = [1   1   1   1
        2   0   0   0
        3   0   0   0];
    
nodal_disp = [1    0   0   0
              2    0   0   0
              3    0   0   0];

% ------------------- SECOND EXAMPLE ---------------------- 

% num_ele = 5/0.1;         % number of elements 
% j = 1 : num_ele;
% conn = [j'   ones(num_ele,6)];          
% supp = [1   1   1   0;  [j(2:end)' zeros(num_ele-1,3)]; num_ele+1 0 1 0];    
% nodal_disp = [j' zeros(num_ele,3); num_ele+1 0 0 0]; 
                    
% ------------------ THIRD EXAMPLE ---------------------- 

% num_ele = 21;
% j = 1:num_ele;
% conn = [j' ones(num_ele,6)];               
% supp = [j(1:end-1)' zeros(num_ele-1,3)];
% supp(1,2:end) = [1 1 1];
% supp(8,2:end) = [1 1 1];
% nodal_disp = [j(1:end-1)' zeros(num_ele-1,3)];

% --------------------- FOURTH EXAMPLE -----------------------

% conn = [1   1   1   1   1   1   1
%         2   1   1   1   1   1   1
%         3   1   1   1   1   1   1
%         4   1   1   1   1   1   1
%         5   1   1   1   1   1   1
%         6   1   1   1   1   1   1
%         7   1   1   1   1   1   1
%         8   1   1   1   1   1   1
%         9   1   1   1   1   1   1
%         10   1   1   1   1   1   1
%         11   1   1   1   1   1   1
%         12   1   1   1   1   1   1];   
%     
% supp =  [1   0   0   0
%          2   0   0   0
%          3   0   0   0
%          4   0   0   0
%          5   0   0   0
%          6   0   0   0
%          7   0   0   0
%          8   0   0   0
%          9   0   0   0
%          10   0   0   0
%          11   0   0   0
%          12   0   0   0
%          13   1   1   1];  
%      
% nodal_disp = [1   0   0   0
%               2   0   0   0
%               3   0   0   0
%               4   0   0   0
%               5   0   0   0
%               6   0   0   0
%               7   0   0   0
%               8   0   0   0
%               9   0   0   0
%              10   0   0   0
%              11   0   0   0
%              12   0   0   0
%              13   0   0   0];  
  
% --------------------- FIFTH EXAMPLE -----------------------
% j = 1:25;
% conn = [j' ones(25,6)];     
% supp =  [[j 26]' zeros(26,3)];
% supp([1,6,11,16,21,26],2:3) = ones(6,2);
% nodal_disp = [[j 26]' zeros(26,3)];     


end