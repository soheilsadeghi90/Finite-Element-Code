function T = trans_mat(Theta)

% produce rotation matrix 

n_ele = length(Theta);

% Initialization
T = zeros(n_ele,6,6);

% v = [1 0 0];

for i = 1:n_ele
    T(i,:,:) = [cos(Theta(i))      sin(Theta(i))      0       0               0           0;
               -sin(Theta(i))      cos(Theta(i))      0       0               0           0;
                0               0               1       0               0           0;
                0               0               0       cos(Theta(i))      sin(Theta(i))  0;
                0               0               0      -sin(Theta(i))      cos(Theta(i))  0;
                0               0               0       0               0           1];
end

end