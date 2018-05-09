function K_glb = glb_stf(be,T,Ke,n_ele,n_nd)

% This function assemble the global stiffness matrix by adding up all
% element matrices in a careful manner, using boolean and rotation
% matrices, according to the given equation.

% Initialization
K_glb = zeros(3*n_nd);

for i = 1:n_ele
    temp = squeeze(be(i,:,:))'*squeeze(T(i,:,:))'*squeeze(Ke(i,:,:))*squeeze(T(i,:,:))*squeeze(be(i,:,:));
    K_glb = K_glb + temp;
end
