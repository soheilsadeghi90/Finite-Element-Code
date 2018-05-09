function K = FE_timo(len_ele,ele_mat,area_ele,areaS_ele,I_ele)

% Generating stiffness matrix for a FE Timoshenko linear shape function
% element

% The matrix is produced as the integral of B*C*B'*jacobian. First, all
% parts are defined parameterically:

syms xi
N = [0.5 * (1 - xi) 0.5 * (1 + xi)];
dif_N = diff(N);

E = ele_mat(1);
G = ele_mat(1)/2/(1+ele_mat(2));
A = area_ele;
As = areaS_ele;
I = I_ele;
L = len_ele;

j_inv = 2/L;

B = [j_inv*dif_N(1)    0                 0               j_inv*dif_N(2)  0                 0;
          0            j_inv*dif_N(1)   -N(1)            0               j_inv*dif_N(2)   -N(2);
          0            0                 j_inv*dif_N(1)  0               0                 j_inv*dif_N(2)];

C = [E*A    0        0;
     0      G*As     0;
     0      0        E*I];

 % Assembling the integrant
integrant = B'*C*B/j_inv;

% Finding the gaussian nodes and weights 
[x,w] = lgwt(2,-1,1);

% Initialization
K = zeros(6,6);

% Integration using Gaussian Quadrature 
for i = 1:length(x)
    K = K + double(subs(integrant,xi,x(i)))*w(i);
end

end

