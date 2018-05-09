function [M_fnl, K_fnl, u0_glb, u_dot0_glb, P_fnl, indx] = rmv_0_stiff(M_fnl, K_fnl, P_fnl, u0_glb, u_dot0_glb,varargin)

% This function is responsible to remove DOFs with zero stiffness from
% mass, stiffness and also intial conditions. It also sends out the indeces
% for the ill DOFs to be used for restoration after integrations. 

indx = fliplr(find(~(diag(K_fnl))));

for i = indx;
    M_fnl(i,:) = [];
    M_fnl(:,i) = [];
    if nargin == 5
        u0_glb(i) = [];
        u_dot0_glb(i) = [];
    end
    P_fnl(i) = [];
end

indx = fliplr(indx);

if nargin == 3
    u0_glb = [];
    u_dot0_glb = [];
end


end

