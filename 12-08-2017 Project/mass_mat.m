function [m_diag] = mass_mat(ele_type, ele_mat, len_ele, area_ele, I_ele, areaS_ele, diag_mtd, idxb, idxa)

% In this function, mass matrix is calculated. The method of calculation
% depends on type of element and diagonalization method which should
% be entered by the program operator. 

E = ele_mat(1);
v = ele_mat(2);
ro = ele_mat(3);
L = len_ele;
I = I_ele;
A = area_ele;
G = E/2/(1+v);
m_ele = zeros(6,6);
m_diag = zeros(6,6);
As = areaS_ele;

syms x
N = [(L-x)/L    x/L];
integ = ro*N'*N*A;
m_ax = int(integ,x,0,L);
me = A*L*ro;

    switch ele_type
         case 1
             m_ele([1,4],[1,4]) = m_ax;
             m_diag = m_ele;

         case 2     % generate mass matrix based on the direct equation for the Timoshenko Beam (Archer 1963). In its current format, it gives very large values!
                beta = E*I/L/L/G/As;
                m_flx = me*[48*beta^2+(42/5)*beta+13/35    (-6*beta^2-11*beta/10-11/210)*L    24*beta^2+18/5*beta+9/70    (6*beta^2+9/10*beta+13/420)*L
                             (-6*beta^2-11*beta/10-11/210)*L   (6*beta^2/5+beta/5+1/105)*L^2   (-6*beta^2-9*beta/10-13/420)*L   (-6*beta^2/5-beta/5-1/140)*L^2
                             24*beta^2+18*beta/5+9/70   (-6*beta^2-9*beta/10-13/420)*L   48*beta^2+42*beta/5+13/35  (6*beta^2+11*beta/10+11/210)*L
                             (6*beta^2+9*beta/10+13/420)*L   (-6*beta^2/5-beta/5-1/140)*L^2   (6*beta^2+11*beta/10+11/210)*L   (6*beta^2/5+beta/5+1/105)*L^2] + ...
                             ro*I/L*[6/5    (6*beta-1/10)*L    6/5  (6*beta-1/10)*L
                             (6*beta-1/10)*L   (48*beta^2+2*beta+2/15)*L^2   (-6*beta+1/10)*L   (24*beta^2-2*beta-1/30)*L^2
                             6/5    (-6*beta+1/10)*L    6/5  (-6*beta+1/10)*L
                             (6*beta-1/10)*L   (24*beta^2-2*beta-1/30)*L^2   (-6*beta+1/10)*L   (48*beta^2+2*beta+2/15)*L^2];
                 m_ele([2,3,5,6],[2,3,5,6]) = m_flx;
                 m_ele([1,4],[1,4]) = m_ax;
                 m_diag = m_ele;                 

         case {3,21}  
              if strcmp(diag_mtd,'cons-SF') == 1    % Mass matrix generation using Shape Function integration. 
                 m_flx = zeros(4,4);
                 syms x
                 N = [1-3*(x^2)/(L^2)+2*(x^3)/(L^3); x-2*(x^2)/L+(x^3)/(L^2); 3*(x^2)/(L^2)-2*(x^3)/(L^3); -(x^2)/L+(x^3)/(L^2)];
                 intg = ro*A*N*N';
                 [x_n,w]=lgwt(2,0,L);
                 for i = 1:length(x_n)
                     m_flx = m_flx + double(subs(intg,x,x_n(i)))*w(i);
                 end
                 m_ax = zeros(2,2);
                 syms x;
                 N = [1-x/L x/L];
                 intg = ro*A*N*N';
                 [x_n,w]=lgwt(2,0,L);
                 for i = 1:length(x_n)
                     m_ax = m_ax + double(subs(intg,x,x_n(i)))*w(i);
                 end
                 m_ele([2,3,5,6],[2,3,5,6]) = m_flx;
                 m_ele([1,4],[1,4]) = m_ax; 
                 m_diag = m_ele;
              else              % Mass matrix generation for the Bernoulli beam using the direct mass matrix fomulation.
                 m_flx = me/420*[156    22*L    54    -13*L
                                 22*L   4*L^2   13*L  -3*L^2
                                 54     13*L    156   -22*L
                                 -13*L  -3*L^2  -22*L 4*L^2];
                 m_ele([2,3,5,6],[2,3,5,6]) = m_flx;
                 m_ele([1,4],[1,4]) = m_ax;
                 m_diag = m_ele;
             end
    end      
% Up to this line, the initial matrix which is consistent has been
% produced. Below, based on the diagonalization method selected by the
% operator, a diagonal matrix will be generated from the consistent mass
% matrix. 
    switch diag_mtd % In case of 'None', consistent mass matrix will be used
        case 'None'
            m_diag = m_ele;
        case 'row summing'      % row summing method for diagonalization
            m_diag = diag(sum(m_ele,2));
        case 'HRZ'              % HRZ method for diagonalization
            x = sum(sum(m_ele([1,4],[1,4])))/(m_ele(1,1)+m_ele(4,4));
            y = sum(sum(m_ele([2,5],[2,5])))/(m_ele(2,2)+m_ele(5,5));
            z = sum(sum(m_ele([3,6],[3,6])))/(m_ele(3,3)+m_ele(6,6));
            m_diag = [m_ele(1,1)*x    0    0    0    0    0
                      0    m_ele(2,2)*y    0    0    0    0
                      0    0    m_ele(3,3)*z    0    0    0
                      0    0    0    m_ele(4,4)*x    0    0
                      0    0    0    0    m_ele(5,5)*y    0
                      0    0    0    0    0    m_ele(6,6)*z];
        case 'ad hoc'          % ad hoc method for diagonalization  
          m_ele = zeros(6,6);
          m_diag = zeros(6,6);  
          switch ele_type
               case 1
                     m_ele = 0.5 * me *[1  0   0   0   0   0
                                        0  0   0   0   0   0
                                        0  0   0   0   0   0
                                        0  0   0   1   0   0
                                        0  0   0   0   0   0
                                        0  0   0   0   0   0];
                 case {2,3,21}
                     m_trans = 0.5 * me * [1   0   0   0
                                           0   1   0   0
                                           0   0   1   0
                                           0   0   0   1];
                     m_rot = me * L^2 / 24 * [1  0
                                              0  1];
                     m_ele([1,2,4,5],[1,2,4,5]) = m_trans;
                     m_ele([3,6],[3,6]) = m_rot;
                     m_diag = m_ele;
          end
    end


% if ele_type ~= 1
%    m_aa = m_diag(idxa,idxa);
%    m_bb = m_diag(idxb,idxb);
%    m_ab = m_diag(idxa,idxb);
%    m_con = m_aa-m_ab*inv(m_bb)*m_ab';
%    m_cond(idxa,idxa) = m_con;
% end
n_a = length(idxa);
n_b = length(idxb);
   
% if ele_type ~= 1
%    m_aa = m_diag(idxa,idxa);
%    m_bb = m_diag(idxb,idxb);
%    m_ab = m_diag(idxa,idxb);
%    T = [eye(n_a,n_a)       zeros(n_a,n_b)
%        -inv(m_bb)*m_ab'    zeros(n_b,n_b)];
%    m_diag = T'*m_diag*T;
% end

end