function [K_fnl, P_fnl] = BC_adj(K_glb,P_glb,supp,nodal_disp)

% Adjust global matrix and load vector with respect to the nodal supports.
% The function zeros all values on the associated col and row if it is
% constrained and overwrite 1 for the associated diagonal value. It also
% modifies load vector by nodal displacements, if exist. 

% Initialization
K_fnl = K_glb;
P_fnl = P_glb;

supp = supp(:);
nodal_disp = nodal_disp(:);

n = length(supp);

for i = 1:n;
    if supp(i) == 1
        K_fnl(i,:) = 0;
        K_fnl(:,i) = 0;
        K_fnl(i,i) = 1;
        P_fnl = P_fnl - K_glb(:,i) * nodal_disp(i);
        P_fnl(i) = nodal_disp(i);
    end
    if supp(i) == 0
        P_fnl = P_fnl - K_glb(:,i) * nodal_disp(i);
    end
    if K_fnl(i,i) == 0
        K_fnl(i,i) = 1;
    end
end

for i = 1:n;
    if supp(i) == 1
        P_fnl(i) = nodal_disp(i);
    end
end


end


        