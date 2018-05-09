function new_mat = addnode(mat,indx)
% This function adds removed nodes in a careful manner, in a way that
% finally the order of nodes will be the same as the beginning order. 
[k,l] = size(mat);
new_mat = zeros(k+length(indx),l);
new_mat(indx,:) = zeros(length(indx),l);
t = setdiff(1:k+length(indx),indx);
new_mat(t,:) = mat;

end