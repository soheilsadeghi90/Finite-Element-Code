function P_glb = glb_ld(P,T,Ge,be,n_ele,n_nd,tag_nd)

% This function assemble the global nodal load vector by adding up all
% equivalent nodal load vector for elements in a careful manner, using 
% boolean and rotation matrices, according to the given equations.

% Initialization
P_glb = zeros(3*n_nd,1);
Ptemp = zeros(3*n_nd,1);

for i = 1:n_ele
    temp = squeeze(be(i,:,:))'*squeeze(T(i,:,:))'*Ge(i,:)';
    P_glb = P_glb + temp;
end

for i = 1 : length(P)
    for j = 1 : n_nd
        if tag_nd(j) == i && j >= i
            Ptemp((j-1)*3+1) = P(i,1);
            Ptemp((j-1)*3+2) = P(i,2);
            Ptemp(j*3) = P(i,3);
        end
    end
end

P_glb = P_glb + Ptemp;

end
