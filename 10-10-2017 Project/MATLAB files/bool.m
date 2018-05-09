function be = bool(st_node,end_node,tag_node,n_ele)

% Produce boolean matrix to map local DOFs to global DOFs for each element.

n_nd = length(tag_node);
be = zeros(n_ele,6,3*n_nd);

for i = 1:n_ele
    for j = 1:n_nd
        if st_node(i) == tag_node(j)
            be(i,1,(j-1)*3+1) = 1;
            be(i,2,(j-1)*3+2) = 1;
            be(i,3,(j-1)*3+3) = 1;
        end
        if end_node(i) == tag_node(j)
            be(i,4,(j-1)*3+1) = 1;
            be(i,5,(j-1)*3+2) = 1;
            be(i,6,(j-1)*3+3) = 1;
        end            
    end
end

end
