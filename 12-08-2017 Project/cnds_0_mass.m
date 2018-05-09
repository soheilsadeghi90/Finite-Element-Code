function [M_fnl, K_fnl, u0_glb, u_dot0_glb, P_fnl,T_k, idxa] = cnds_0_mass(M_fnl, K_fnl, P_fnl, u0_glb, u_dot0_glb, varargin)

% % This function performs the matrix condensation on all parts of the 
% equation of motion, such as Mass matrix, Stiffness matrix, loading and 
% intial conditions vectors, using the transformation matrix method. Having
% the transformation matrix, decondensation process very is straight forward.
% NOTE: all of the matrices mentioned above will be condensated by the
% transformation matrix produced by stiffness matrix. However, the ill DOFs
% that get condensated are those that have zero mass. 

idxb = find(~(diag(M_fnl)));
idxa = find((diag(M_fnl)));
n_a = length(idxa);
n_b = length(idxb);

k_aa = K_fnl(idxa,idxa);
k_bb = K_fnl(idxb,idxb);
k_ab = K_fnl(idxa,idxb);
T_k = [eye(n_a,n_a)       zeros(n_a,n_b)        % Stiffness transformation matrix
    -inv(k_bb)*k_ab'    zeros(n_b,n_b)];

K_fnl = T_k'*K_fnl*T_k;
P_fnl = T_k'*P_fnl;
M_fnl = T_k'*M_fnl*T_k;

if nargin == 5
    u0_glb = T_k'*u0_glb;
    u_dot0_glb = T_k'*u_dot0_glb;
end

% m_aa = M_fnl(idxa,idxa);
% m_bb = M_fnl(idxb,idxb);
% m_ab = M_fnl(idxa,idxb);
% T_m = [eye(n_a,n_a)       zeros(n_a,n_b)
%     -inv(m_bb)*m_ab'    zeros(n_b,n_b)];
% M_fnl = T_m'*M_fnl*T_m;

K_fnl = real(K_fnl+K_fnl')/2;

if nargin == 3
    u0_glb = [];
    u_dot0_glb = [];
end



end

