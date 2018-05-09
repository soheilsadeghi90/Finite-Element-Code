function [u0_glb,u_dot0_glb] = glb_ic(u0,u_dot0,T,be,n_ele,n_nd)

% This function assemble the global nodal load vector by adding up all
% equivalent nodal load vector for elements in a careful manner, using 
% boolean and rotation matrices, according to the given equations.

u0_glb = u0';
u_dot0_glb = u_dot0';

u0_glb = u0_glb(:);
u_dot0_glb = u_dot0_glb(:);

% for i = 1:n_ele
%     temp_u = squeeze(be(i,:,:))'*squeeze(T(i,:,:))'*u0(i,:);
%     u0_glb = u0_glb + temp;
%     temp_udot = squeeze(be(i,:,:))'*squeeze(T(i,:,:))'*u_dot0(i,:);
%     u_dot0_glb = u_dot0_glb + temp;
% end
% 
% for i = 1 : length(P)
%     for j = 1 : n_nd
%         if tag_nd(j) == i && j >= i
%             Ptemp((j-1)*3+1) = P(i,1);
%             Ptemp((j-1)*3+2) = P(i,2);
%             Ptemp(j*3) = P(i,3);
%         end
%     end
% end

% P_glb = P_glb + Ptemp;

end
