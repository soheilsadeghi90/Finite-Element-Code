function [ge,ge_a,ge_b,idxb,idxa] = Eq_nd_ld(len_ele,DL_ele,conn)

% produce equivalent nodal load vectors for each element

% Initialization
ge = zeros(6,1);

L = len_ele;
% Distributed axial load
ge_1 = [DL_ele(1)*L/2  0	0	DL_ele(1)*L/2   0	0]';
% Distributed transverse load
ge_2 = [0	DL_ele(2)*L/2 DL_ele(2)*L^2/12   0	DL_ele(2)*L/2  -DL_ele(2)*L^2/12]';
% Distributed moment load
ge_3 = [0	0	DL_ele(3)*L/2	0	0	DL_ele(3)*L/2]';

ge = ge_1 +  ge_2 + ge_3;

idxa = [];
idxb = [];

% Finding disconnected DOFs for each element and save its label for
% condensation (Discontinuity resolution in 'stf_mat.m')

for j = 1:6
    if conn(j) == 0
       idxb = [idxb j];
    end
    if conn(j) == 1
       idxa = [idxa j];
    end        
end

ge_b = ge(idxb);
ge_a = ge(idxa);

end
