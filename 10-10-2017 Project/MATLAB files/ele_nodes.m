function [tag_ele,tag_nd,st_node,end_node,Len,ele_prop_ID_n,DL,ele_typ_n,MAT_ID_n,conn_n,Theta] = ...
    ele_nodes(x, y,ele_type,MAT_ID,conn,ele_prop_ID)

% This function is a two sided function. First is an INPUT function. So,
% the operator needs to enter inputs. For the five verifying examples,
% inputs are prepared. Uncomment it if needed.

% Maximum size of mesh needed (the unit is the same as units used for other inputs)
mesh_max = 10;                % Maximum size of the Mesh

% ----------------- First Example ------------------------

% ele_num = 50;                % number of elements to model the nonprismatic beam
% tag_nd = 1:length(x);        
% st_node = [1:ele_num];       % starting nodes for elements
% end_node = [2:ele_num + 1];  % ending nodes for elements 
% DL = zeros(ele_num,3);       % Distributed Load for elements - row per element (parallel, transverse, moment)

% -------------------- SECOND EXAMPLE ---------------------

% tag_nd = 1:length(x);        
% st_node =  [1   2];          % starting nodes for elements
% end_node = [2   3];          % ending nodes for elements   
% DL = [0 0   0;0 0   0];      % Distributed Load for elements - row per element (parallel, transverse, moment)

% -------------------- THIRD EXAMPLE ---------------------

% tag_nd = 1:length(x);
% n_nd = length(x);
% st_node = [1    3   5   2   4   6   3   5   7];   % starting nodes for elements
% end_node = [3   5   7   4   6   8   4   6   8];   % ending nodes for elements  
% 
% dist_l = 6000; % [N/cm]
% DL = [0 0   0;
%     0   0   0
%     0   0   0
%     0   0   0
%     0   0   0
%     0   0   0
%     0   0   0
%     0   0   0
%     0   -dist_l   0];       % Distributed Load for elements - row per element (parallel, transverse, moment)


% -------------------- FOURTH EXAMPLE ---------------------

% st_node =  [1   2   3   4   5   6   7   8   8   7   7   7   6];   % starting nodes for elements
% end_node = [2   3   4   5   6   7   8   1   2   2   3   4   4];   % ending nodes for elements  
% 
% DL = [   0   0   0
%          0   0   0
%          0   0   0
%          0   0   0
%          0   0   0
%          0   0   0
%          0   0   0
%          0   0   0
%          0   0   0
%          0   0   0
%          0   0   0
%          0   0   0
%          0   0   0      ];    % Distributed Load for elements - row per element (parallel, transverse, moment)
     
% -------------------- FIFTH EXAMPLE ---------------------

st_node =  [1   2   3   4];         % starting nodes for elements
end_node = [2   3   4   5];         % ending nodes for elements    

DL = [   0   0       0
         0   0       0
         0   -1000   0
         0   0       0      ];      % Distributed Load for elements - row per element (parallel, transverse, moment)     

     
% -------------------- SIXTH EXAMPLE ---------------------

% st_node =  [1   2];         % starting nodes for elements
% end_node = [2   3];         % ending nodes for elements    
% 
% DL = [   0   0       0
%          0   0       0
%          0   0       0
%          0   0       0      ];      % Distributed Load for elements - row per element (parallel, transverse, moment)     
     
%% print element locations

% as a part of input that is printed

disp('Elements Locations:');
fprintf('\n')
for i = 1:length(st_node)
    disp(['element ',num2str(i),' connects node ',num2str(st_node(i)),' to node ',num2str(end_node(i))]);
end
fprintf('\n\n\n')


%% Find length and angle of elements

n_ele = length(st_node);

for i = 1:n_ele
    Len(i) = sqrt((x(st_node(i))-x(end_node(i)))^2+(y(st_node(i))-y(end_node(i)))^2);
    Theta(i) = atan((y(end_node(i)) - y(st_node(i)))/(x(end_node(i)) - x(st_node(i))));
end

tag_ele = 1:n_ele;

%% MESHING SUB-SECTION

i = 1;

% The meshing is done in this while loop. It checks it an element is longer
% than the max size of mesh. If so, cut it to pieces and modify all other
% parameters correspondingly. For example, if an element is splitted to 10
% pieces, nine nodes are placed between start and end nodes of that
% element. Correspondingly, distributed loads, connectivities and other
% props are duplicated/updated in a careful manner. 

while i <= n_ele
% for i = 1:n_ele
    if ele_type(i) == 21
        if Len(i) > mesh_max
            n_add_nd = ceil(Len(i)/mesh_max) - 1;
            mesh_ele_ln = Len(i)/(n_add_nd + 1);
            j = st_node(i)*1000 + end_node(i)*100 + 1:st_node(i)*1000 + end_node(i)*100 + n_add_nd; % multiplier is very sensitive! (100 works well)
            st_node = [st_node(1:i) j st_node(i+1:end)];
            end_node = [end_node(1:i-1) j end_node(i:end)];
            Len = [Len(1:i-1) mesh_ele_ln*ones(1,n_add_nd + 1) Len(i+1:end)];
            Theta = [Theta(1:i-1) Theta(i)*ones(1,n_add_nd + 1) Theta(i+1:end)];
            ele_type = [ele_type(1:i-1) ele_type(i)*ones(1,n_add_nd + 1) ele_type(i+1:end)];
            ele_prop_ID = [ele_prop_ID(1:i-1) ele_prop_ID(i)*ones(1,n_add_nd + 1) ele_prop_ID(i+1:end)];
            DL = [DL(1:i-1,:); repmat(DL(i,:),n_add_nd + 1,1);DL(i+1:end,:)];
            tag_ele = [tag_ele(1:i-1) j j(end)+1 tag_ele(i+1:end)];
            n_ele = length(st_node);
            MAT_ID = [MAT_ID(1:i-1) MAT_ID(i)*ones(1,n_add_nd + 1) MAT_ID(i+1:end)];
            conn = [tag_ele(1:i-1)' conn(1:i-1,2:end); tag_ele(i) conn(i,2:4) 1 1 1; ...
                tag_ele(i+1:i+n_add_nd-1)' repmat([1 1 1 1 1 1],length(tag_ele(i+1:i+n_add_nd-1)),1); ...
                tag_ele(i+n_add_nd)' 1 1 1 conn(i,5:7); tag_ele(i+n_add_nd+1:end)' conn(i+1:end,2:end)];
        end
    end
    i = i + 1;
end

ele_typ_n = ele_type;
MAT_ID_n = MAT_ID;
conn_n = conn;
tag_nd = union(st_node,end_node);      
ele_prop_ID_n = ele_prop_ID;

end

