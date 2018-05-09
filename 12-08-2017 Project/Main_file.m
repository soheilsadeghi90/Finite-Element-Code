%% Main file
clear; clc;

%%
%       ----------------- I N P U T -------------------%

% Commentary: In this section, the operator has to enter all required
% inputs in three sub-functions: 'nodal_coordinate.m', 'ele_nodes.m' and
% 'conn.m'. More descriptions on each inputs are provided on those files. 

% General Analysis Inputs
anlz_typ = 3;                   % type of analysis which is going to be done (1 for static, 2 for modal, 3 for time history)
n_mode = 4;                     % number of modes needed to be comupted or considered for modal and time history analysis
diag_mtd = 'cons-DIR';            % mass matrix diagonalization method ('cons-DIR', 'cons-SF', 'HRZ', 'row summing' or 'ad hoc')
dmp_typ = 'Rayleigh';           % damping type ('None', 'Rayleigh' or 'Modal')

% if Rayleigh damping is chosen, enter alfa and beta
f1 = 10*6.28;   % Min frequency for Rayleigh calculations
f2 = 120*6.28;  % Max frequency for Rayleigh calculations
w1 = f1/2/pi;
w2 = f2/2/pi;
dmp1_R = .02;
dmp2_R = .02;
alfa_R = 2*w1*w2*(dmp1_R*w2-dmp2_R*w1)/(w2^2-w1^2);
beta_R = 2*(dmp2_R*w2-dmp1_R*w1)/(w2^2-w1^2);
% alfa_R = 0.5;     % Just for example 3
% beta_R = 0.5;     % Just for example 3

% if modal damping is chosen, modal damping should be given here
dmp_mdl = 0.05;

% --------------- Loading time history defenition ---------------
% Example one - no load
del_t = 0.001;
t = 0:del_t:0.5;
ld_mltp = zeros(1,length(t));

% Exmaple two - Ramp loading
% del_t = 0.0001;
% t = 0:del_t:0.1;
% ld_interval = 0.01;
% ld_mltp = [0.0:del_t/ld_interval:1 zeros(1,length(t)-length(0.0:del_t/ld_interval:1))];

% Example three - Ramp step loading
% del_t = 0.001;
% t = 0:del_t:1.0;
% ld_interval_1 = 0.1;
% ld_interval_2 = 0.2;
% ld_mltp = [0:1/(ld_interval_1/del_t+1):1 ones(1,ceil((ld_interval_2-ld_interval_1)/del_t))  zeros(1,length(t)-(ld_interval_2/del_t+1)-1)];

% Example Four - Sine loading
% del_t = 0.005;
% t = 0:del_t:5.0;
% ld_interval = 0.5;
% ld_mltp = [sin((0:del_t:ld_interval).*2.*pi./0.5) zeros(1,length(t)-length(0:del_t:ld_interval))];

% Example Five - Multiple loading/continuous beam
% del_t = 0.0005;
% t = 0:del_t:3.0;
% ld_interval = 0.5;         
% ld_mltp = [sin((0:del_t:ld_interval).*2.*pi./0.5) zeros(1,length(t)-length(0:del_t:ld_interval))];
    
% nodal_coordinate calls all nodal coordinates and loads
[tag_nd_orig, x_nd, y_nd, P, u0, u_dot0] = nodal_coordinate();       % tag_nd_orig = node tags as the operator model (It may change if meshing happens)
ele_num = [];                                                        % an initial value needed for non-prismatic sections

% There are also three more inputs that have to be entered in this file
% (main.m). Those are material_ID, element_type and element_prop_ID which
% call different properties from the built-in libraries. Notes on libraries
% are provided in each of them. 
% For each validation examples, corresponding inputs are prepared and the
% operator just needs to uncomment them:

% ------------ FIRST Example ---------------

material_ID = [1,1];                          % Enter IDs using labels in 'material.m' (material library)
element_type = [2,2];                         % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)   
element_prop_ID = [4,4];                      % Enter prop_ID labels from sec_prop.m (mechanical props library)

% ------------ SECOND Example ---------------

% ele_num = 5/0.1;
% material_ID = 2*ones(1,ele_num);              % Enter IDs using labels in 'material.m' (material library)
% element_type = 3*ones(1,ele_num);               % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)    
% element_prop_ID = 4*ones(1,ele_num);          % Enter prop_ID labels from sec_prop.m (mechanical props library)

% ------------ THIRD Example --------------

% material_ID = 3*ones(1,21);                          % Enter IDs using labels in 'material.m' (material library)
% element_type = 3*ones(1,21);                       % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)   
% element_prop_ID = 1*ones(1,21);                      % Enter prop_ID labels from sec_prop.m (mechanical props library)

% ------------ FOURTH Example ---------------

% material_ID = 4*ones(1,12);                   % Enter IDs using labels in 'material.m' (material library)                      
% element_type = 2*ones(1,12);                  % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)   
% element_prop_ID = 2*ones(1,12);               % Enter prop_ID labels from sec_prop.m (mechanical props library)

% ------------ FIFTH Example --------------

% material_ID = 4*ones(1,25);                   % Enter IDs using labels in 'material.m' (material library)                      
% element_type = 2*ones(1,25);                  % Enter element types (1: BAR, 2: TIMOSHENKO, 3: BERNOULLI, 21: FEM-TIMO)   
% element_prop_ID = 3*ones(1,25);               % Enter prop_ID labels from sec_prop.m (mechanical props library)

% -------------------------------------------

% The last part is the same for all examples

% function conn.m calls all connectivities, nodal supports and nodal displacement
[conn_o,supp,nodal_disp] = conn();  

% Newmark Parameters. If anlz_typ is 3, parameters are necessary to be entered
acc_asump = 1;      % acceleration assumption for Newmark (1 for constant, 2 for linear, 3 for manual parameters)
switch acc_asump
    case 1
        gamma_N = 0.5;
        beta_N = 0.25;
    case 2
        gamma_N = 0.5;
        beta_N = 1/6;        
        if del_t > (1/f2)*(dmp_mdl*(gamma_N-0.5)+sqrt(0.5*gamma_N-beta_N+dmp_mdl^2*(gamma_N-0.5)^2))/(gamma_N-beta_N)
            display ('del_t is not small enough to stabilize the Newmark method');
            return
        end
    case 3
        gamma_N = 0.6;      % If Manual, enter gamma here
        beta_N = 0.3025;    % If Manual, enter beta here 
        if gamma_N > 2*beta_N || gamma_N < 0.5
            if del_t > (1/f2)*(dmp_mdl*(gamma_N-0.5)+sqrt(0.5*gamma_N-beta_N+dmp_mdl^2*(gamma_N-0.5)^2))/(gamma_N-beta_N)
                display ('del_t is not small enough to stabilize the Newmark method');
                return
            end
        end
end

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

[tag_ele,tag_nd,ele_st_nd,ele_end_nd,len_ele,ele_prop_ID,DL_ele,ele_type,MAT_ID,conn,Theta,u0,u_dot0] ...
    = ele_nodes(x_nd, y_nd, element_type,material_ID,conn_o,element_prop_ID,u0,u_dot0); 

%%
%       ----------------- A N A L Y S I S -------------------%

% store number of nodes and elements (it is done after meshing)
n_nd = length(tag_nd);
n_ele = length(tag_ele);


%% 
%       ------------------- C H E C K S --------------------%

% check the length of ld_mltp wrt time vector

if anlz_typ == 3
    if length(ld_mltp) ~= length(t)
        display ('loading time history and time vector are not the same size, so inconsistent!');
        return
    end
end

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
Me = zeros(n_ele,6,6);
T_stf = zeros(n_ele,6,6);
T_mass = zeros(n_ele,6,6);

if anlz_typ == 1 || anlz_typ == 2

    % LOOP OVER ELEMENTS. produce stiffness and mass matrices and nodal loads
    for i = 1:n_ele
        % load
        [ge,ge_a,ge_b,idxb,idxa] = Eq_nd_ld(len_ele(i),DL_ele(i,:),conn(i,2:7)); 
        if anlz_typ == 2
            [Me(i,:,:)] = mass_mat(ele_type(i), ele_mat(i,:), len_ele(i), area_ele(i), I_ele(i), areaS_ele(i), diag_mtd,idxb,idxa);
            [Me(i,:,:),T_mass(i,:,:)] = condens(Me(i,:,:),idxb,idxa);
        end
        % stiffness matrix
        [Ke(i,:,:),Ge(i,:)] = stf_mat(len_ele(i),ele_mat(i,:),ele_type(i),area_ele(i),I_ele(i),areaS_ele(i),idxb,idxa,ge_a,ge_b,ge);
        [Ke(i,:,:),T_stf(i,:,:)] = condens(Ke(i,:,:),idxb,idxa);
    end
    
    % transformation matrix
    T = trans_mat(Theta);
    
    % Boolean Matrix
    be = bool(ele_st_nd,ele_end_nd,tag_nd,n_ele);
    
    % Global Stiffness Matrix
    K_glb = glb_stf(be,T,Ke,n_ele,n_nd);
    if anlz_typ == 2
        M_glb = glb_stf(be,T,Me,n_ele,n_nd);
    end
    
    % Global Loading Vector
    P_glb = glb_ld(P,T,Ge,be,n_ele,n_nd,tag_nd);
    
    % Impose Boundary Conditions
    if anlz_typ == 1 
        M_glb = zeros(size(K_glb));
        [M_fnl, K_fnl, P_fnl, u0_fnl, u_dot0_fnl,indx_BC] = BC_adj(K_glb,P_glb,supp(:,2:4)',nodal_disp(:,2:4)', M_glb);
    % Removing 0_stiffness dofs from K, M, u0, v0, and P
        [M_fnl, K_fnl, u0_fnl, u_dot0_fnl, P_fnl, indx_stf0] = rmv_0_stiff(M_fnl, K_fnl, P_fnl);  
    else
        [M_fnl, K_fnl, P_fnl, u0_fnl, u_dot0_fnl,indx_BC] = BC_adj(K_glb,P_glb,supp(:,2:4)',nodal_disp(:,2:4)', M_glb);
    % Removing 0_stiffness dofs from K, M, u0, v0, and P
        [M_fnl, K_fnl, u0_fnl, u_dot0_fnl, P_fnl, indx_stf0] = rmv_0_stiff(M_fnl, K_fnl, P_fnl);   
    % Condensate 0_mass dofs from K,M, u0, v0, and P
        [M_fnl, K_fnl, u0_fnl, u_dot0_fnl, P_fnl, T_k, indxa_m0] = cnds_0_mass(M_fnl, K_fnl, P_fnl);        
    end

end

if anlz_typ == 3

    % LOOP OVER ELEMENTS. produce stiffness and mass matrices and nodal loads
    for i = 1:n_ele
        % load
        [ge,ge_a,ge_b,idxb,idxa] = Eq_nd_ld(len_ele(i),DL_ele(i,:),conn(i,2:7)); 
        % Mass Matrix
        [Me(i,:,:)] = mass_mat(ele_type(i), ele_mat(i,:), len_ele(i), area_ele(i), I_ele(i), areaS_ele(i), diag_mtd,idxb,idxa);
        % Stiffness Matrix
        [Ke(i,:,:),Ge(i,:)] = stf_mat(len_ele(i),ele_mat(i,:),ele_type(i),area_ele(i),I_ele(i),areaS_ele(i),idxb,idxa,ge_a,ge_b,ge);
    end
    
    % transformation matrix
    T = trans_mat(Theta);
    
    % Boolean Matrix
    be = bool(ele_st_nd,ele_end_nd,tag_nd,n_ele);
    
    % Global Stiffness Matrix
    K_glb = glb_stf(be,T,Ke,n_ele,n_nd);
    M_glb = glb_stf(be,T,Me,n_ele,n_nd);
    
    
    % Global Loading Vector
    P_glb = glb_ld(P,T,Ge,be,n_ele,n_nd,tag_nd);
    
    % Global Initial Conditions
    [u0_glb,u_dot0_glb] = glb_ic(u0,u_dot0,T,be,n_ele,n_nd);    
    
    % Impose Boundary Conditions
    [M_fnl, K_fnl, P_fnl, u0_fnl, u_dot0_fnl, indx_BC] = BC_adj(K_glb,P_glb,supp(:,2:4)',nodal_disp(:,2:4)', M_glb, u0_glb, u_dot0_glb);
    
    % Removing 0_stiffness dofs from K, M, u0, v0, and P
    [M_fnl, K_fnl, u0_fnl, u_dot0_fnl, P_fnl, indx_stf0] = rmv_0_stiff(M_fnl, K_fnl, P_fnl, u0_fnl, u_dot0_fnl);
    
    % Condensate 0_mass dofs from K,M, u0, v0, and P
    [M_fnl, K_fnl, u0_fnl, u_dot0_fnl, P_fnl, T_k, indxa_m0] = cnds_0_mass(M_fnl, K_fnl, P_fnl, u0_fnl, u_dot0_fnl);
    
    % Damping Matrix
    switch dmp_typ
        case 'None'
            C_fnl = zeros(size(K_fnl));
        case 'Rayleigh'
            C_fnl = alfa_R * M_fnl + beta_R * K_fnl;
        case 'Modal'
            vec = dmp_mdl*ones(1,size(K_fnl,1));
            C_fnl = diag(vec);
    end
    
 
    % ---------------------- NEWMARK MATHOD ----------------------
    
    Ntsteps = length(ld_mltp);
    % response time history initialization
    u_t = zeros(length(u0_fnl),Ntsteps);
    v_t = zeros(length(u0_fnl),Ntsteps);
    a_t = zeros(length(u0_fnl),Ntsteps);
    u_t(:,1) = u0_fnl;
    v_t(:,1) = u_dot0_fnl;
    a_t(:,1) = M_fnl\(P_fnl-C_fnl*u_dot0_fnl-K_fnl*u0_fnl);
    K_eff = 1/(beta_N * del_t^2) * M_fnl + gamma_N / (beta_N * del_t) * C_fnl + K_fnl;
    
    for i = 1:Ntsteps-1
        F_hat = P_fnl*ld_mltp(i) + M_fnl*((1/beta_N/del_t/del_t)*u_t(:,i) + (1/beta_N/del_t)*v_t(:,i) + (1/2/beta_N-1)*a_t(:,i)) ...
            + C_fnl*((gamma_N/beta_N/del_t)*u_t(:,i) + (gamma_N/beta_N - 1)*v_t(:,i) + del_t*(gamma_N/2/beta_N - 1)*a_t(:,i));
        u_t(:,i+1) = K_eff\F_hat;
        a_t(:,i+1) = (1/beta_N/del_t/del_t)*(u_t(:,i+1)-u_t(:,i)-del_t*v_t(:,i))-(1/2/beta_N - 1)*a_t(:,i);
        v_t(:,i+1) = gamma_N/beta_N/del_t*(u_t(:,i+1)-u_t(:,i))-(gamma_N/beta_N - 1)*v_t(:,i)-del_t*(gamma_N/2/beta_N - 1)*a_t(:,i);
    end
  
% restore 0-mass,0-stifness, and BC's
    u_t = T_k\u_t;
    v_t = T_k\v_t;
    a_t = T_k\a_t;
    stf_res_u = addnode(u_t,indx_stf0);
    BC_res_u = addnode(stf_res_u, indx_BC);
    v_t = T_k\v_t;
    stf_res_v = addnode(v_t,indx_stf0);
    BC_res_v = addnode(stf_res_v, indx_BC);
    a_t = T_k\a_t;
    stf_res_a = addnode(a_t,indx_stf0);
    BC_res_a = addnode(stf_res_a, indx_BC);    
    
end

%%       
%       ------------------- P O S T  P R O C E S S --------------------%

switch anlz_typ
    case 1
        % solve governing equation for displacements
        U = K_fnl\P_fnl;
        U_res_mod = addnode(U,indx_stf0);
        BC_res_U = addnode(U_res_mod, indx_BC);        
        % find nodal forces
        F_c = K_glb*BC_res_U - P_glb;
    case 2
        % solve eigenvalue analysis (MODAL ANALYSIS)
        [modeShape,freq] = eigs(K_fnl,M_fnl,n_mode,'sm');
        freq = sqrt(diag(freq))/(2*pi);
        norm_modal_m = diag( modeShape' * M_fnl * modeShape).^(1/2);
        norm_modes = modeShape * diag(1./norm_modal_m);
        % restore DOF's
        modeShape = T_k\modeShape;
        stf_res_mod = addnode(modeShape,indx_stf0);
        BC_res_mod = addnode(stf_res_mod, indx_BC);
        Fcmodes = K_fnl * norm_modes - M_fnl * norm_modes * diag(freq)./2./pi;
end

if anlz_typ == 2
    % Mode Shapes plot
    prompt = 'Which mode shape do you want to plot? ';
    mode_id_to_plot = input(prompt);
    prompt2 = 'displacement scale factor? ';
    Mode = BC_res_mod(:,end+1-mode_id_to_plot)';
    scl_fct = input(prompt2);
    Mode = Mode.*scl_fct;
    x_new = x_nd + [Mode(1:3:end)];
    y_new = y_nd + [Mode(2:3:end)];
    for i = 1:length(ele_st_nd)
        plot([x_nd(ele_st_nd(i)), x_nd(ele_end_nd(i))],[y_nd(ele_st_nd(i)),y_nd(ele_end_nd(i))],'b'); hold on;
    end
    hold on;
    for i = 1:length(ele_st_nd)
        plot([x_new(ele_st_nd(i)), x_new(ele_end_nd(i))],[y_new(ele_st_nd(i)),y_new(ele_end_nd(i))],'r'); hold on;
    end
    frq = fliplr(freq');
    title(['Mode ',num2str(mode_id_to_plot), ', scale factor = ', num2str(scl_fct), ', frequency = ', num2str(frq(mode_id_to_plot))])
end


% displaying outputs (displacements and forces at each node)
% NOTE that secondary nodes produced by meshing are not considered for display
switch anlz_typ
    case 1
        disp('------------- DISPLACEMENTS AT NODES -------------');
        fprintf('\n')
        for i = 1:length(tag_nd_orig);
            disp(['X displacement at node ' num2str(i),' = ', num2str(BC_res_U((i-1)*3+1))]);
            disp(['Y displacement at node ' num2str(i),' = ', num2str(BC_res_U((i-1)*3+2))]);
            disp(['Z Rotation at node     ' num2str(i),' = ', num2str(BC_res_U((i-1)*3+3))]);
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
    case 2
        fprintf('\n')
        disp('------------- MODAL FREQUENCIES -------------');
        frq = fliplr(freq');
        for i = 1:n_mode;
            disp(['Mode ' num2str(i),' frequency = ', num2str(frq(i))]);
            fprintf('\n')
        end  
    case 3
        prompt = 'Which node displacement do you want to plot? ';
        d1 = input(prompt);
        prompt2 = 'Which DOF? (1,2 or 3) ';
        d2 = input(prompt2);
        plot(t,BC_res_u((d1-1)*3+d2,:));
end
        








