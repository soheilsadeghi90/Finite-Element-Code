%% Main file
clear; clc;

%%
%       ----------------- I N P U T -------------------%

% Commentary: In this section, the operator has to enter all required
% inputs in three sub-functions: 'nodal_coordinate.m', 'ele_nodes.m' and
% 'conn.m'. More descriptions on each inputs are provided on those files. 

% nodal_coordinate calls all nodal coordinates and loads
[tag_nd_orig, x_nd, y_nd, P] = nodal_coordinate();       % tag_nd_orig = node tags as the operator model (It may change if meshing happens)
ele_num = [];                                            % an initial values needed for non-prismatic sections

% There are also three more inputs that have to be entered in this file
% (main.m). Those are material_ID, element_type and element_prop_ID which
% call different properties from the built-in libraries. Notes on libraries
% are provided in each of them. 
% For each validation examples, corresponding inputs are prepared and the
% operator just needs to uncomment them:

% ------------ First Example ---------------

% ele_num = 50;
% material_ID = 2*ones(1,ele_num);              % Enter IDs using labels in 'material.m' (material library)
% element_type = ones(1,ele_num);               % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)    
% element_prop_ID = 3*ones(1,ele_num);          % Enter prop_ID labels from sec_prop.m (mechanical props library)

% ------------ Second Example --------------

% material_ID = [1,1];                          % Enter IDs using labels in 'material.m' (material library)
% element_type = [21,21];                       % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)   
% element_prop_ID = [4,4];                      % Enter prop_ID labels from sec_prop.m (mechanical props library)

% ------------ Third Example ---------------

% material_ID = [3,3,3,3,3,3,3,3,3];            % Enter IDs using labels in 'material.m' (material library)                      
% element_type = [2,2,2,2,2,2,2,2,2];           % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)   
% element_prop_ID = [1,1,1,1,1,1,2,2,2];        % Enter prop_ID labels from sec_prop.m (mechanical props library)

% ------------ Fourth Example --------------

% material_ID = ones(13,1);                     % Enter IDs using labels in 'material.m' (material library)
% element_type = 1*ones(13,1);                  % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)  
% element_prop_ID = 5 * ones(13,1);             % Enter prop_ID labels from sec_prop.m (mechanical props library)

% ------------ Fifth Example ---------------

material_ID = [3,3,3,3];                        % Enter IDs using labels in 'material.m' (material library)
element_type = [2,2,2,2];                       % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)     
element_prop_ID = [2,2,2,2];                    % Enter prop_ID labels from sec_prop.m (mechanical props library)

% ------------ Sixth Example ---------------

% material_ID = [3,3];                        % Enter IDs using labels in 'material.m' (material library)
% element_type = [21,21];                       % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)     
% element_prop_ID = [2,2];                    % Enter prop_ID labels from sec_prop.m (mechanical props library)

% -------------------------------------------

% The last part is the same for all examples

% function conn.m calls all connectivities, nodal supports and nodal displacement
[conn_o,supp,nodal_disp] = conn();  

% print inputs up to this point:

disp('------------------ INPUTS -----------------');
fprintf('\n\n')
disp('Nodal Coordinates:');
fprintf('\n')
for i = 1:length(tag_nd_orig)
    disp(['node ',num2str(i),' : [',num2str(x_nd(i)),',',num2str(y_nd(i)),']']);
end
fprintf('\n\n')

disp('Connectivity:');
fprintf('\n')
disp('element i : [x1  y1  rot1  x2  y2  rot2]; 1 if connected and 0 if not')
fprintf('\n')
for i = 1:length(material_ID)
    disp(['element ',num2str(i),' : [',num2str(conn_o(i,2:7)),']']);
end
fprintf('\n\n')

disp('Constraints:');
fprintf('\n')
disp('node i contrains: [x y rot]; 1 if yes and 0 if no')
fprintf('\n')
for i = 1:length(tag_nd_orig)
    disp(['node ',num2str(i),' constraints: ',num2str(supp(i,2:end))]);
end
fprintf('\n\n')

disp('Nodal Displacements:');
fprintf('\n')
disp('node i displacements: [dx dy drot]; enter values if any')
fprintf('\n')
for i = 1:length(tag_nd_orig)
    disp(['node ',num2str(i),' displacements: ',num2str(nodal_disp(i,2:end))]);
end
fprintf('\n\n')

disp('Material IDs: (material definitions can be found in material.m file)');
fprintf('\n')
disp('element i material ID : ID number')
fprintf('\n')
for i = 1:length(material_ID)
    disp(['element ',num2str(i),' material ID : ',num2str(material_ID(i))]);
end
fprintf('\n\n')

disp('Element types: (1 if BAR / 2 if TIMOSHENKO / 3 if BERNOULLI / 21 if FE TIMO LINEAR SF)');
fprintf('\n')
disp('element i type : type label')
fprintf('\n')
for i = 1:length(material_ID)
    disp(['element ',num2str(i),' type : ',num2str(element_type(i))]);
end
fprintf('\n\n')

disp('Element Section Property ID');
fprintf('\n')
disp('element i section ID : (section definitions can be found in sec_prop.m file)')
fprintf('\n')
for i = 1:length(material_ID)
    disp(['element ',num2str(i),' section ID : ',num2str(element_type(i))]);
end
fprintf('\n\n')

% function ele_nodes has different parts inside. 1) load start and end
% nodes for each element, name elements and find lengths and angles for
% each of them. 2) If Finite Elements is chosen, it meshes large elements
% automatically. Notet that for meshing, a max mesh length should be
% entered in the ele_nodes.m file.

[tag_ele,tag_nd,ele_st_nd,ele_end_nd,len_ele,ele_prop_ID,DL_ele,ele_type,MAT_ID,conn,Theta] ...
    = ele_nodes(x_nd, y_nd, element_type,material_ID,conn_o,element_prop_ID); 

%%
%       ----------------- A N A L Y S I S -------------------%

% store number of nodes and elements (it is done after meshing)
n_nd = length(tag_nd);
n_ele = length(tag_ele);


%% 
%       ------------------- C H E C K S --------------------%

% check if conn.m input data are entered sufficiently.
if size(conn,1) ~= n_ele
    display ('Connectivity data are not compatible with the number of elements!');
    return
end

if size(supp,1) ~= tag_nd_orig
    display ('Support data are not compatible with the number of nodes!');
    return
end
if size(nodal_disp,1) ~= tag_nd_orig
    display ('Nodal displacement data are not compatible with the number of elements!');
    return
end

% checking loads with respect to the element type. No transverse or moment
% loads on BAR element.

if sum(abs(DL_ele(:,2))) ~= 0 || sum(abs(DL_ele(:,3))) ~= 0 
    for i = 1:n_ele
        if ele_type(i) == 1
            display ('Transverse load or moment is applied on Bar Element!');
            return
        end
    end
end

% check number of material IDs and element types. should be same as the number of elements

if length(MAT_ID) ~= n_ele
    display('Material IDs are not the same as the number of elements!');
    return
end
if length(ele_type) ~= n_ele
    display('Element types entered are not the same as the number of elements!');
    return
end

% call material props from library [material.m]
[ele_mat] = material(MAT_ID);

% call mechanical props per section from library [sec_prop.m]
[I_ele, area_ele, areaS_ele] = sec_prop(ele_prop_ID, ele_num);  

% Initialization for local stiffness matrix and eq. nodal load vec
Ke = zeros(n_ele,6,6);
Ge = zeros(n_ele,6);

% LOOP OVER ELEMENTS. produce stiffness matrices and nodal loads
for i = 1:n_ele
    % load
    [ge,ge_a,ge_b,idxb,idxa] = Eq_nd_ld(len_ele(i),DL_ele(i,:),conn(i,2:7)); 
    % stiffness matrix
    [Ke(i,:,:),Ge(i,:)] = stf_mat(len_ele(i),ele_mat(i,:),ele_type(i),area_ele(i),I_ele(i),areaS_ele(i),idxb,idxa,ge_a,ge_b,ge);
end

% transformation matrix
T = trans_mat(Theta);

% Boolean Matrix
be = bool(ele_st_nd,ele_end_nd,tag_nd,n_ele);
    
% Global Stiffness Matrix
K_glb = glb_stf(be,T,Ke,n_ele,n_nd);

% Global Loading Vector
P_glb = glb_ld(P,T,Ge,be,n_ele,n_nd,tag_nd);

% Impose Boundary Conditions
[K_fnl, P_fnl] = BC_adj(K_glb,P_glb,supp(:,2:4)',nodal_disp(:,2:4)');

%%       
%       ------------------- P O S T  P R O C E S S --------------------%

% solve governing equation for displacements
U = K_fnl\P_fnl;

% find nodal forces
F_c = K_glb*U - P_glb;

% displaying outputs (displacements and forces at each node)
% NOTE that secondary nodes produced by meshing are not considered for display

disp('------------- DISPLACEMENTS AT NODES -------------');
fprintf('\n')
for i = 1:length(tag_nd_orig);
    disp(['X displacement at node ' num2str(i),' = ', num2str(U((i-1)*3+1))]);
    disp(['Y displacement at node ' num2str(i),' = ', num2str(U((i-1)*3+2))]);
    disp(['Z Rotation at node     ' num2str(i),' = ', num2str(U((i-1)*3+3))]);
    fprintf('\n')
end

disp('------------- FORCES AT NODES -------------');
fprintf('\n')
for i = 1:length(tag_nd_orig);
    disp(['X force at node ' num2str(i),' = ', num2str(F_c((i-1)*3+1))]);
    disp(['Y force at node ' num2str(i),' = ', num2str(F_c((i-1)*3+2))]);
    disp(['moment at node  ' num2str(i),' = ', num2str(F_c((i-1)*3+3))]);
    fprintf('\n')
end








