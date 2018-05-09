function [M_fnl, K_fnl, P_fnl, u0_fnl, u_dot0_fnl,indx_BC] = BC_adj(K_glb,P_glb,supp,nodal_disp,M_glb, u0_glb, u_dot0_glb,varargin) 

% Adjust global matrix and load vector with respect to the nodal supports.
% The function zeros all values on the associated col and row if it is
% constrained and overwrite 1 for the associated diagonal value. It also
% modifies load vector by nodal displacements, if exist. 
% Initialization

if nargin == 7
    M_fnl = M_glb;
    u0_fnl = u0_glb;
    u_dot0_fnl = u_dot0_glb;
elseif nargin == 6
    return
elseif nargin == 5
    M_fnl = M_glb;
    u0_fnl = [];
    u_dot0_fnl = [];    
elseif nargin == 4
    M_fnl = [];
    u0_fnl = [];
    u_dot0_fnl = [];    
end

supp = supp(:);
nodal_disp = nodal_disp(:);
n = length(supp);
 
% for i = 1:n;
%     if supp(i) == 1
%         K_fnl(i,:) = 0;
%         K_fnl(:,i) = 0;
%         K_fnl(i,i) = 1;
%         if nargin == 5
%             M_fnl(i,:) = 0;
%             M_fnl(:,i) = 0;
%             M_fnl(i,i) = 1;
%         end
%         P_fnl = P_fnl - K_glb(:,i) * nodal_disp(i);
%         P_fnl(i) = nodal_disp(i);
%     end
%     if supp(i) == 0
%         P_fnl = P_fnl - K_glb(:,i) * nodal_disp(i);
%     end
%     if K_fnl(i,i) == 0
%         K_fnl(i,i) = 1;
%     end
%     if nargin == 5 && M_fnl(i,i) == 0
%         M_fnl(i,i) = 1;
%     end
% end
% 
% for i = 1:n;
%     if supp(i) == 1
%         P_fnl(i) = nodal_disp(i);
%     end
% end

indx = [];
indx_BC = [];

for i = 1:n;
    if supp(i) == 0
        indx = [indx i];
    else
        indx_BC = [indx_BC i];
    end
end

K_fnl = K_glb(indx,indx);
M_fnl = M_glb(indx,indx);
P_fnl = P_glb(indx);
if nargin == 7
    u0_fnl = u0_fnl(indx);
    u_dot0_fnl = u_dot0_fnl(indx);
end

        