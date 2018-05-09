function [K,Ge] = stf_mat(len_ele,ele_mat,ele_type,area_ele,I_ele,areaS_ele,idxb,idxa,ge_a,ge_b,ge)

K = zeros(6,6);

% ele_type = 1 : Bar
% ele_type = 2 : Timoshenko - Direct Matrix
% ele_type = 3 : Bernoulli - Direct Matrix
% ele_type = 21: FE Linear Shape Function [Timoshenko]
   
E = ele_mat(1);
v = ele_mat(2);
L = len_ele;
A = area_ele;
G = E/2/(1+v);
I = I_ele;
As = areaS_ele;

switch ele_type
     case 1
            K = E*A/L*[ 1   0   0   -1  0   0
                        0   0   0   0   0   0
                        0   0   0   0   0   0
                        -1  0   0   1   0   0
                        0   0   0   0   0   0
                        0   0   0   0   0   0];
     case 2
            phi = 12*E*I/L/L/G/As;
            K = 1/(1+phi)*[E*A/L*(1+phi),  0,           0,            -E*A/L*(1+phi),  0,             0
                           0,               12*E*I/L^3,  6*E*I/L^2,     0,               -12*E*I/L^3,   6*E*I/L^2
                           0,               6*E*I/L^2,  (4+phi)*E*I/L,  0,               -6*E*I/L^2,    (2-phi)*E*I/L
                           -E*A/L*(1+phi),  0,           0,             E*A/L*(1+phi),  0,             0
                           0,              -12*E*I/L^3, -6*E*I/L^2,     0,              12*E*I/L^3,    -6*E*I/L^2
                           0,               6*E*I/L^2,   (2-phi)*E*I/L, 0,              -6*E*I/L^2,    (4+phi)*E*I/L];
     case 3
            K =            [E*A/L       0              0            -E*A/L      0               0;
                            0           12*E*I/L^3	   6*E*I/L^2     0         -12*E*I/L^3,     6*E*I/L^2;
                            0           6*E*I/L^2      4*E*I/L       0         -6*E*I/L^2	    2*E*I/L;
                           -E*A/L       0              0             E*A/L      0               0;
                            0          -12*E*I/L^3	  -6*E*I/L^2     0          12*E*I/L^3     -6*E*I/L^2;
                            0           6*E*I/L^2	   2*E*I/L       0         -6*E*I/L^2       4*E*I/L]; 
    case 21
           % FE_timo.m function generates FE linear shape function for
           % timoshenko beam.         
           K = FE_timo(len_ele,ele_mat,area_ele,areaS_ele,I_ele);
end

%% Discontinuity - stiffness and nodal load vector modifications

% discontinuities are resolved here by equating zero connections by ones.
% NOTE that this part is not implemented on BAR elements. 

n_a = length(idxa);
n_b = length(idxb);

% if ele_type ~= 1
%    k_aa = K(idxa,idxa);
%    k_bb = K(idxb,idxb);
%    k_ab = K(idxa,idxb);
%    k_con = k_aa-k_ab*inv(k_bb)*k_ab';
%    if size(k_ab*inv(k_bb)*ge_b,2) == 0
%        g_con = ge_a;
%    else
%       g_con = ge_a - k_ab*inv(k_bb)*ge_b;
%    end
%    K(idxa,idxa) = k_con;
%    ge(idxa) = g_con;
%    Ge = ge;
% end

% if ele_type ~= 1
%    k_aa = K(idxa,idxa);
%    k_bb = K(idxb,idxb);
%    k_ab = K(idxa,idxb);
%    T = [eye(n_a,n_a)       zeros(n_a,n_b)
%        -inv(k_bb)*k_ab'    zeros(n_b,n_b)];
%    K = T'*K*T;
%    Ge = T'*ge;
% end

Ge = ge;



            
            
