function [M,T] = condens(M,idxb,idxa)

% This function performs the matrix condensation using the transformation
% matrix method. Having the transformation matrix, decondensation process
% is straight forward. 

M = squeeze(M);

n_a = length(idxa);
n_b = length(idxb);

% if ele_type ~= 1
m_aa = M(idxa,idxa);
m_bb = M(idxb,idxb);
m_ab = M(idxa,idxb);
T = [eye(n_a,n_a)       zeros(n_a,n_b)
    -inv(m_bb)*m_ab'    zeros(n_b,n_b)];
M = T'*M*T;
% end

end