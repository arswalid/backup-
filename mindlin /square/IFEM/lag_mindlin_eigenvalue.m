
clear
clc
clf
column=1;
%  25 35 45 55
for n_elems_x_o = [190]-1 %odd number  bigger or equal to 3

format long 
E1=5.e6;
E2=5.e6;
% these value are for a compression shock at 0.

reduced =0;

nu=  0.3;
k = 5/6;
Lx= 10;
Ly= 10;
thickness = Lx*0.1; 

rho=1;

mode = 4;

C_s1 = ((E1*thickness*k)/(2*(1+nu))) * [1 0 ; 0 1 ];

C_b1 = ((E1*thickness^3)/(12*(1-nu^2))) * [1 nu 0;nu 1 0;0 0 (1-nu)/2 ];


C_s2 = ((E2*thickness*k)/(2*(1+nu))) * [1 0 ; 0 1 ];

C_b2 = ((E2*thickness^3)/(12*(1-nu^2))) * [1 nu 0;nu 1 0;0 0 (1-nu)/2 ];
%domain size


n_elems_y_o = n_elems_x_o;   
n_nodes_x_o = n_elems_x_o+1;
n_nodes_y_o = n_elems_y_o+1;
n_elems_o   = n_elems_x_o*n_elems_y_o;
n_nodes_o   = n_nodes_x_o*n_nodes_y_o;
n_dofs_o    = n_nodes_o;
h_x       = Lx / n_elems_x_o;
h_y       = Ly / n_elems_y_o;

n_elems_x = n_elems_x_o+1; %odd number   %interface
n_elems_y = n_elems_y_o;
n_nodes_x = n_elems_x+1+1;  %2 interface
n_nodes_y = n_elems_y+1;
n_elems   = n_elems_x*n_elems_y;
n_nodes   = n_nodes_x*n_nodes_y;
n_dofs    = n_nodes;%mult 4
%
% global indices are stored using the following rule
%  u   values for all nodes at y=0
%  followed by u values at all nodes at y=h_y, etc.
%

% indices for boundary nodes
bottom_edge_indices   =   (1:n_nodes_x);
top_edge_indices      =   (n_nodes-n_nodes_x+1:n_nodes);
left_edge_indices     =   (1:n_nodes_x:n_nodes-n_nodes_x+1);
right_edge_indices    =   (n_nodes_x:n_nodes_x:n_nodes);

% nodal locations to plot numerical solution)
x_vals = [(0:h_x:(h_x*((n_elems_x_o/2)-0.5)))';((n_elems_x_o/2)+0.5)*h_x;(h_x*((n_elems_x_o/2)-0.5));...
    ((n_elems_x_o/2)+0.5)*h_x;(((n_elems_x_o/2)+1.5)*h_x:h_x:Lx)'] ;
y_vals = (0:h_y:Ly)';
x_loc  = zeros(n_nodes,2);
n_elem_nodes = 4;
for i=1:n_nodes_y
    x_loc((i-1)*n_nodes_x+1:i*n_nodes_x,1) = x_vals;
    x_loc((i-1)*n_nodes_x+1:i*n_nodes_x,2) = y_vals(i);
end


% four point Gauss-quadrature rule for 2D domain integration
qp_loc_domain  = ...
    [-1/sqrt(3)  1/sqrt(3)  -1/sqrt(3)  1/sqrt(3); ...  % xi location
    -1/sqrt(3)  -1/sqrt(3)  1/sqrt(3)  1/sqrt(3)];     % eta location
qp_wgt_domain  = [1 1 1 1];
% two point Gauss-quadrature rule for bottom edge of elem
qp_loc_bottom  = ...
    [-1/sqrt(3)  1/sqrt(3); ...  % xi location
    -1         -1 ];     % eta location
qp_loc_bottom1  = ...
    [-0.7887  -0.2113; ...  % xi location
    -1         -1 ];     % eta location
qp_loc_bottom2  = ...
    [0.2113     0.7887; ...  % xi location
    -1         -1 ];     % eta location
% two point Gauss-quadrature rule for top edge of elem
qp_loc_top     = ...
    [-1/sqrt(3)  1/sqrt(3); ...  % xi location
    1          1 ];     % eta location
qp_loc_top1     = ...
    [-0.7887  -0.2113; ...  % xi location
    1          1 ];     % eta location
qp_loc_top2     = ...
    [0.2113     0.7887; ...  % xi location
    1          1 ];     % eta location
% two point Gauss-quadrature rule for left edge of elem
qp_loc_left    = ...
    [      -1           -1; ...  % xi location
    -1/sqrt(3)  1/sqrt(3)];      % eta location
% two point Gauss-quadrature rule for right edge of elem
qp_loc_right   = ...
    [       1            1; ...  % xi location
    -1/sqrt(3)  1/sqrt(3)];      % eta location
% four point Gauss-quadrature rule for 2D domain integration
qp_loc_domainhalf1  = ...
    [-0.7887  -0.2113  -0.2113  -0.7887; ...  % xi location
    -1/sqrt(3)  -1/sqrt(3)  1/sqrt(3)  1/sqrt(3)];     % eta location

qp_wgt_domain1  = [1/2 1/2 1/2 1/2];

qp_loc_domainhalf2  = ...
    [0.2113     0.7887    0.7887 0.2113  ; ...  % xi location
    -1/sqrt(3)  -1/sqrt(3)  1/sqrt(3)  1/sqrt(3)];     % eta location

qp_wgt_domain2  = [1/2 1/2 1/2 1/2];


% two point Gauss-quadrature rule for disc
qp_loc_disc   = ...
    [       0            0; ...  % xi location
    -1/sqrt(3)  1/sqrt(3)];      % eta location


qp_wgt_boundary1  = [1 1];
qp_wgt_boundary2  = [1/2 1/2];
n_qp_domain  = 4;
n_qp_boundary= 2;



% lag_bc1          = zeros(n_dofs+1, n_dofs+1);
% lag_bc2          = zeros(n_dofs+1, n_dofs+1);
% lag_bc3          = zeros(n_dofs+1, n_dofs+1);
% lag_bc4          = zeros(n_dofs+1, n_dofs+1);
lag     = sparse(n_dofs*3, 3*n_dofs);

X_vec   = zeros(3*n_dofs, 1);
dX_vec  = zeros(3*n_dofs, 1);

elem1=1;

lag_ind_elem1=1;
lag_ind_elem2=1;

C_b = C_b1;
C_s = C_s1;
for i_elem_x=1:n_elems_x+1
    for i_elem_y=1:n_elems_y
        
        i_elem                  = (i_elem_y-1)*n_elems_x+i_elem_x;
        indices_x               = [i_elem_x i_elem_x+1];
        indices_y               = [i_elem_y i_elem_y+1];
        indices                 = ...
            [(indices_y(1)-1)*n_nodes_x+indices_x ...
            (indices_y(2)-1)*n_nodes_x+fliplr(indices_x)];
         w_vec                   = X_vec(indices);   % initial solution estimate
        tetax                   = X_vec(indices+n_nodes);
        tetay                   = X_vec(indices+2*n_nodes); 
        x_vec                   = x_loc(indices,1);
        y_vec                   = x_loc(indices,2);
        if i_elem_x == n_elems_x/2
            
            % zero all matrices before evaluation
            n_elem_dofs     = max(size(indices));
            jac_e           = zeros(n_elem_dofs*3,n_elem_dofs*3);
       
            res_e           = zeros(n_elem_dofs*3,n_elem_dofs*3);
            
            
            for i_qp = 1:n_qp_domain
                
                xi                      = qp_loc_domainhalf1(1,i_qp);           % element quadrature point location
                eta                     = qp_loc_domainhalf1(2,i_qp);           % element quadrature point location
                
                % solution and derivative quantities at quadrature
                % point
[w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);

            A = [ zeros(1,4) zeros(1,4) -Bxmat;...
                  zeros(1,4) Bymat   zeros(1,4);...
                  zeros(1,4) Bxmat      -Bymat];
             
             
            D = [Bmat zeros(1,4) zeros(1,4)];
                
                % calculate the stiff  and mass matrix
            res_e                   = res_e + qp_wgt_domain1(i_qp) * J_det *(...
                       ...                   % -dphi/dx du/dx
                      ...                   % -dphi/dy du/dy
                + D'  * D);                              % +phi f 
            
            jac_e                   = jac_e + qp_wgt_domain1(i_qp) * J_det *(...
                A' * C_b1 * A       ...                   
              );       
            end
            
            
        if reduced == 1
        
                                for i_qp = 1
                
                xi                      = 0;           % element quadrature point location
                eta                     = 0;           % element quadrature point location
                
                % solution and derivative quantities at quadrature
                % point
[w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);

  
              
            B = [Bxmat zeros(1,4) -Bmat;...
                 Bymat  Bmat zeros(1,4)];

            jac_e                   = jac_e + 4 * J_det *(...
                     ...                   
               B' * C_s1 * B );       
                        end
        
        else
        
        for i_qp = 1:n_qp_domain
            
            xi                      = qp_loc_domainhalf1(1,i_qp);           % element quadrature point location
            eta                     = qp_loc_domainhalf1(2,i_qp);           % element quadrature point location
                
                % solution and derivative quantities at quadrature
                % point
[w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);

  
              
            B = [Bxmat zeros(1,4) -Bmat;...
                 Bymat  Bmat zeros(1,4)];

            jac_e                   = jac_e + qp_wgt_domain1(i_qp) * J_det *(...
                     ...                   
               B' * C_s1 * B );       
                                        end
        end
            
                    
for i=1:4
    ind((i-1)*4+1:i*4,1) = indices;
    ind((i-1)*4+1:i*4,2) = indices(i);
end

ind1=ind+n_dofs;
ind2=ind+2*n_dofs;


z1=zeros(1,16);
z2=zeros(1,16);
z3=zeros(1,16);
z4=zeros(1,16);
z5=zeros(1,16);
z6=zeros(1,16);
z7=zeros(1,16);
z8=zeros(1,16);
z9=zeros(1,16);

o=1;
p=1;
q=1;
col = 0;
for j =1:12
       if j == 5 || j == 9
            o=1;
            p=1;
            q=1;
        end 
    for i =1:12
        if i>=1 && i<=4 && j>=1 && j<=4
            z1(o+col) = jac_e(i,j);
            o=o+1;
        elseif i>=5 && i<=8 && j>=1 && j<=4
            z2(p+col) = jac_e(i,j);
            p=p+1; 
        elseif i>=9 && i<=12 && j>=1 && j<=4
            z3(q+col) = jac_e(i,j);
            q=q+1;
        elseif i>=1 && i<=4 && j>=5 && j<=8
            z4(o+col) = jac_e(i,j);
            o=o+1;
        elseif i>=5 && i<=8 && j>=5 && j<=8
            z5(p+col) = jac_e(i,j);
            p=p+1;
        elseif i>=9 && i<=12 && j>=5 && j<=8
            z6(q+col) = jac_e(i,j);
            q=q+1;
        elseif i>=1 && i<=4 && j>=9 && j<=12
            z7(o+col) = jac_e(i,j);
            o=o+1;
        elseif i>=5 && i<=8 && j>=9 && j<=12
            z8(p+col) = jac_e(i,j);
            p=p+1;
        elseif i>=9 && i<=12 && j>=9 && j<=12
            z9(q+col) = jac_e(i,j);
            q=q+1;
        end
 
    end

end
clear o p q 


vec_i_el = [ind(:,1);ind1(:,1);ind2(:,1);ind(:,1);ind(:,1);ind1(:,1);ind1(:,1);ind2(:,1);ind2(:,1)];
vec_j_el = [ind(:,2);ind1(:,2);ind2(:,2);ind1(:,2);ind2(:,2);ind(:,2);ind2(:,2);ind(:,2);ind1(:,2)];
z_el     = [z1';z5';z9';z4';z7';z2';z8';z3';z6'];

if elem1 == 1 
vec_i_stiff = [ vec_i_el];
vec_j_stiff = [ vec_j_el];
z_stiff=[z_el];
else
vec_i_stiff = [vec_i_stiff ; vec_i_el];
vec_j_stiff = [vec_j_stiff ; vec_j_el];
z_stiff=[z_stiff;z_el];
end
clear z_el vec_i vec_j

   o=1;
for j = 1:4
    for i=1:4
z(o) = res_e(i,j);
      o=o+1;
    end
end
   clear o 
   
   vec_i_el = [ind(:,1);ind1(:,1);ind2(:,1)];
   vec_j_el = [ind(:,2);ind1(:,2);ind2(:,2)];
   
   z_el  = [(rho*thickness)*z'  ;  ((rho*thickness^3)/12) *z' ;  ((rho*thickness^3)/12) *z'  ] ;
   if elem1 == 1 
   vec_i_mass = [ vec_i_el];
   vec_j_mass = [ vec_j_el];
   z_mass=[z_el];  
   else
   vec_i_mass = [vec_i_mass ; vec_i_el];
   vec_j_mass = [vec_j_mass ; vec_j_el];
   z_mass=[z_mass;z_el];
   end
   clear z vec_i vec_j

            %lagrange
            lag1           = zeros(n_elem_dofs,1);
            lag2           = zeros(n_elem_dofs,1);
            lag3           = zeros(n_elem_dofs,1);
            for i_qp = 1:n_qp_boundary
                
                xi                      = qp_loc_disc(1,i_qp);           % element quadrature point location
                eta                     = qp_loc_disc(2,i_qp);           % element quadrature point location
                
                
                % solution and derivative quantities at quadrature
                % point
                % right edge uses mode = 2
              [w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 2);
                
                % calculate the residual and Jacobians
                lag1           =  lag1  + qp_wgt_boundary1(i_qp) * J_det *(Bmat');
                lag2           =  lag2  + qp_wgt_boundary1(i_qp) * J_det *(Bmat');
                lag3           =  lag3  + qp_wgt_boundary1(i_qp) * J_det *(-Bmat');
            end
            
            lag             (          indices,3*n_dofs   +lag_ind_elem1  )=lag1;
            lag             (n_dofs   +indices,3*n_dofs   +lag_ind_elem1+1)=lag2;
            lag             (2*n_dofs +indices,3*n_dofs   +lag_ind_elem1+2)=lag3;
            
            lag             (3*n_dofs    +lag_ind_elem1  ,            indices)=lag1';
            lag             (3*n_dofs    +lag_ind_elem1+1,   n_dofs  +indices)=lag2';
            lag             (3*n_dofs    +lag_ind_elem1+2,   2*n_dofs+indices)=lag3';
            
            lag_ind_elem1=lag_ind_elem1+3;
        elseif i_elem_x == n_elems_x/2+2
            
            C_b = C_b2;
            C_s = C_s2;
            
            n_elem_dofs     = max(size(indices));
            jac_e           = zeros(n_elem_dofs*3,n_elem_dofs*3);
            res_e           = zeros(n_elem_dofs*3,n_elem_dofs*3);
            for i_qp = 1:n_qp_domain
                
                xi                      = qp_loc_domainhalf2(1,i_qp);           % element quadrature point location
                eta                     = qp_loc_domainhalf2(2,i_qp);           % element quadrature point location
                
                         % solution and derivative quantities at quadrature
                % point
[w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);

            A = [ zeros(1,4) zeros(1,4) -Bxmat;...
                  zeros(1,4) Bymat   zeros(1,4);...
                  zeros(1,4) Bxmat      -Bymat];
             
            D = [Bmat zeros(1,4) zeros(1,4)];
                
                % calculate the stiff  and mass matrix
            res_e                   = res_e + qp_wgt_domain2(i_qp) * J_det *(...
                       ...                   % -dphi/dx du/dx
                      ...                   % -dphi/dy du/dy
                + D'  * D);                              % +phi f 
            
            jac_e                   = jac_e + qp_wgt_domain2(i_qp) * J_det *(...
                A' * C_b * A       ...                   
               );       
            end
        if reduced == 1
        
                                for i_qp = 1
                
                xi                      = 0;           % element quadrature point location
                eta                     = 0;           % element quadrature point location
                
                % solution and derivative quantities at quadrature
                % point
[w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);

  
              
            B = [Bxmat zeros(1,4) -Bmat;...
                 Bymat  Bmat zeros(1,4)];

            jac_e                   = jac_e + 4 * J_det *(...
                     ...                   
               B' * C_s1 * B );       
                        end
        
        else
        
        for i_qp = 1:n_qp_domain
            
            xi                      = qp_loc_domainhalf2(1,i_qp);           % element quadrature point location
            eta                     = qp_loc_domainhalf2(2,i_qp);           % element quadrature point location
                
                % solution and derivative quantities at quadrature
                % point
[w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);

  
              
            B = [Bxmat zeros(1,4) -Bmat;...
                 Bymat  Bmat zeros(1,4)];

            jac_e                   = jac_e + qp_wgt_domain2(i_qp) * J_det *(...
                     ...                   
               B' * C_s1 * B );       
                                        end
        end
for i=1:4
    ind((i-1)*4+1:i*4,1) = indices;
    ind((i-1)*4+1:i*4,2) = indices(i);
end

ind1=ind+n_dofs;
ind2=ind+2*n_dofs;


z1=zeros(1,16);
z2=zeros(1,16);
z3=zeros(1,16);
z4=zeros(1,16);
z5=zeros(1,16);
z6=zeros(1,16);
z7=zeros(1,16);
z8=zeros(1,16);
z9=zeros(1,16);

o=1;
p=1;
q=1;
col = 0;
for j =1:12
       if j == 5 || j == 9
            o=1;
            p=1;
            q=1;
        end 
    for i =1:12
        if i>=1 && i<=4 && j>=1 && j<=4
            z1(o+col) = jac_e(i,j);
            o=o+1;
        elseif i>=5 && i<=8 && j>=1 && j<=4
            z2(p+col) = jac_e(i,j);
            p=p+1; 
        elseif i>=9 && i<=12 && j>=1 && j<=4
            z3(q+col) = jac_e(i,j);
            q=q+1;
        elseif i>=1 && i<=4 && j>=5 && j<=8
            z4(o+col) = jac_e(i,j);
            o=o+1;
        elseif i>=5 && i<=8 && j>=5 && j<=8
            z5(p+col) = jac_e(i,j);
            p=p+1;
        elseif i>=9 && i<=12 && j>=5 && j<=8
            z6(q+col) = jac_e(i,j);
            q=q+1;
        elseif i>=1 && i<=4 && j>=9 && j<=12
            z7(o+col) = jac_e(i,j);
            o=o+1;
        elseif i>=5 && i<=8 && j>=9 && j<=12
            z8(p+col) = jac_e(i,j);
            p=p+1;
        elseif i>=9 && i<=12 && j>=9 && j<=12
            z9(q+col) = jac_e(i,j);
            q=q+1;
        end
 
    end

end
clear o p q 


vec_i_el = [ind(:,1);ind1(:,1);ind2(:,1);ind(:,1);ind(:,1);ind1(:,1);ind1(:,1);ind2(:,1);ind2(:,1)];
vec_j_el = [ind(:,2);ind1(:,2);ind2(:,2);ind1(:,2);ind2(:,2);ind(:,2);ind2(:,2);ind(:,2);ind1(:,2)];
z_el     = [z1';z5';z9';z4';z7';z2';z8';z3';z6'];

if elem1 == 1 
 vec_i_stiff = [ vec_i_el];
vec_j_stiff = [ vec_j_el];
z_stiff=[z_el];
else
vec_i_stiff = [vec_i_stiff ; vec_i_el];
vec_j_stiff = [vec_j_stiff ; vec_j_el];
z_stiff=[z_stiff;z_el];
end
clear z_el vec_i vec_j

   o=1;
for j = 1:4
    for i=1:4
z(o) = res_e(i,j);
      o=o+1;
    end
end
   clear o 
   
   vec_i_el = [ind(:,1);ind1(:,1);ind2(:,1)];
   vec_j_el = [ind(:,2);ind1(:,2);ind2(:,2)];
   
   z_el  = [(rho*thickness)*z'  ;  ((rho*thickness^3)/12) *z' ;  ((rho*thickness^3)/12) *z'  ] ;
   if elem1 == 1 
   vec_i_mass = [ vec_i_el];
   vec_j_mass = [ vec_j_el];
   z_mass=[z_el];  
   else
   vec_i_mass = [vec_i_mass ; vec_i_el];
   vec_j_mass = [vec_j_mass ; vec_j_el];
   z_mass=[z_mass;z_el];
   end
   clear z vec_i vec_j
            
            %lagrange
            lag1           = zeros(n_elem_dofs,1);
            lag2           = zeros(n_elem_dofs,1);
            lag3           = zeros(n_elem_dofs,1);
            
            for i_qp = 1:n_qp_boundary
                
                xi                      = qp_loc_disc(1,i_qp);           % element quadrature point location
                eta                     = qp_loc_disc(2,i_qp);           % element quadrature point location
                
                
                % solution and derivative quantities at quadrature
                % point
                % right edge uses mode = 2
              [w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 2);
                
                % calculate the residual and Jacobians
                lag1           =  lag1  + qp_wgt_boundary1(i_qp) * J_det *(Bmat');
                lag2           =  lag2  + qp_wgt_boundary1(i_qp) * J_det *(Bmat');
                lag3           =  lag3  + qp_wgt_boundary1(i_qp) * J_det *(-Bmat');
            end
            
            
            lag(         indices             ,3*n_dofs   +lag_ind_elem2   )= -lag1;
            lag(n_dofs+  indices             ,3*n_dofs   +lag_ind_elem2 +1)= -lag2;
            lag(2*n_dofs+indices             ,3*n_dofs   +lag_ind_elem2 +2)= -lag3;
            
            lag(3*n_dofs+lag_ind_elem2                 ,  indices)=         -lag1';
            lag(3*n_dofs+lag_ind_elem2+1               ,  indices+  n_dofs)=-lag2';
            lag(3*n_dofs+lag_ind_elem2+2               ,  indices+2*n_dofs)=-lag3';
            
            
         lag_ind_elem2=lag_ind_elem2+3;
            
        elseif i_elem_x == n_elems_x/2+1
            
        else
            
  % zero all matrices before evaluation
        n_elem_dofs     = max(size(indices));
        jac_e           = zeros(n_elem_dofs*3, n_elem_dofs*3);
        res_e           = zeros(n_elem_dofs*3,n_elem_dofs*3);
        
        for i_qp = 1:n_qp_domain
            
            xi                      = qp_loc_domain(1,i_qp);           % element quadrature point location
            eta                     = qp_loc_domain(2,i_qp);           % element quadrature point location
            
            % solution and derivative quantities at quadrature
            % point
[w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);
            
            % calculate the residual and Jacobians
            
            A = [ zeros(1,4) zeros(1,4) -Bxmat;...
                  zeros(1,4) Bymat   zeros(1,4);...
                  zeros(1,4) Bxmat      -Bymat];

             
            D = [Bmat zeros(1,4) zeros(1,4)];
            
            
            res_e                   = res_e + qp_wgt_domain(i_qp) * J_det *(...
                       ...                   % -dphi/dx du/dx
                      ...                   % -dphi/dy du/dy
                + D'  * D);                              % +phi f 
            
            jac_e                   = jac_e + qp_wgt_domain(i_qp) * J_det *(...
                A' * C_b * A       ...                   
                );                         
        end
        
        if reduced == 1
        
                                for i_qp = 1
                
                xi                      = 0;           % element quadrature point location
                eta                     = 0;           % element quadrature point location
                
                % solution and derivative quantities at quadrature
                % point
[w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);

  
              
            B = [Bxmat zeros(1,4) -Bmat;...
                 Bymat  Bmat zeros(1,4)];

            jac_e                   = jac_e + 4 * J_det *(...
                     ...                   
               B' * C_s1 * B );       
                        end
        
        else
        
        for i_qp = 1:n_qp_domain
            
            xi                      = qp_loc_domain(1,i_qp);           % element quadrature point location
            eta                     = qp_loc_domain(2,i_qp);           % element quadrature point location
                
                % solution and derivative quantities at quadrature
                % point
[w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);

  
              
            B = [Bxmat zeros(1,4) -Bmat;...
                 Bymat  Bmat zeros(1,4)];

            jac_e                   = jac_e + qp_wgt_domain(i_qp) * J_det *(...
                     ...                   
               B' * C_s1 * B );       
                                        end
        end
        
for i=1:4
    ind((i-1)*4+1:i*4,1) = indices;
    ind((i-1)*4+1:i*4,2) = indices(i);
end

ind1=ind+n_dofs;
ind2=ind+2*n_dofs;


z1=zeros(1,16);
z2=zeros(1,16);
z3=zeros(1,16);
z4=zeros(1,16);
z5=zeros(1,16);
z6=zeros(1,16);
z7=zeros(1,16);
z8=zeros(1,16);
z9=zeros(1,16);

o=1;
p=1;
q=1;
col = 0;
for j =1:12
       if j == 5 || j == 9
            o=1;
            p=1;
            q=1;
        end 
    for i =1:12
        if i>=1 && i<=4 && j>=1 && j<=4
            z1(o+col) = jac_e(i,j);
            o=o+1;
        elseif i>=5 && i<=8 && j>=1 && j<=4
            z2(p+col) = jac_e(i,j);
            p=p+1; 
        elseif i>=9 && i<=12 && j>=1 && j<=4
            z3(q+col) = jac_e(i,j);
            q=q+1;
        elseif i>=1 && i<=4 && j>=5 && j<=8
            z4(o+col) = jac_e(i,j);
            o=o+1;
        elseif i>=5 && i<=8 && j>=5 && j<=8
            z5(p+col) = jac_e(i,j);
            p=p+1;
        elseif i>=9 && i<=12 && j>=5 && j<=8
            z6(q+col) = jac_e(i,j);
            q=q+1;
        elseif i>=1 && i<=4 && j>=9 && j<=12
            z7(o+col) = jac_e(i,j);
            o=o+1;
        elseif i>=5 && i<=8 && j>=9 && j<=12
            z8(p+col) = jac_e(i,j);
            p=p+1;
        elseif i>=9 && i<=12 && j>=9 && j<=12
            z9(q+col) = jac_e(i,j);
            q=q+1;
        end
 
    end

end
clear o p q 


vec_i_el = [ind(:,1);ind1(:,1);ind2(:,1);ind(:,1);ind(:,1);ind1(:,1);ind1(:,1);ind2(:,1);ind2(:,1)];
vec_j_el = [ind(:,2);ind1(:,2);ind2(:,2);ind1(:,2);ind2(:,2);ind(:,2);ind2(:,2);ind(:,2);ind1(:,2)];
z_el     = [z1';z5';z9';z4';z7';z2';z8';z3';z6'];

if elem1 == 1 
 vec_i_stiff = [ vec_i_el];
vec_j_stiff = [ vec_j_el];
z_stiff=[z_el];
else
vec_i_stiff = [vec_i_stiff ; vec_i_el];
vec_j_stiff = [vec_j_stiff ; vec_j_el];
z_stiff=[z_stiff;z_el];
end
clear z_el vec_i vec_j

   o=1;
for j = 1:4
    for i=1:4
z(o) = res_e(i,j);
      o=o+1;
    end
end
   clear o 
   
   vec_i_el = [ind(:,1);ind1(:,1);ind2(:,1)];
   vec_j_el = [ind(:,2);ind1(:,2);ind2(:,2)];
   
   z_el  = [(rho*thickness)*z'  ;  ((rho*thickness^3)/12) *z' ;  ((rho*thickness^3)/12) *z'  ] ;
   if elem1 == 1 
   vec_i_mass = [ vec_i_el];
   vec_j_mass = [ vec_j_el];
   z_mass=[z_el];  
   else
   vec_i_mass = [vec_i_mass ; vec_i_el];
   vec_j_mass = [vec_j_mass ; vec_j_el];
   z_mass=[z_mass;z_el];
   end
   clear z vec_i vec_j
   
    elem1=elem1+1;
   
        end

    end
end

jac_mat = sparse(vec_i_stiff,vec_j_stiff,z_stiff,3*n_dofs,3*n_dofs);
res_vec = sparse(vec_i_mass,vec_j_mass,z_mass,3*n_dofs,3*n_dofs);


sizeofmatrices=3*n_dofs+lag_ind_elem1-1;
jac_mat(sizeofmatrices,sizeofmatrices)=0;
res_vec(sizeofmatrices,sizeofmatrices)=0;
jac_mat=jac_mat+lag;


jac_mat = jac_mat + speye(sizeofmatrices,sizeofmatrices)*1.e-15;

% res_vec = res_vec + speye(sizeofmatrices,sizeofmatrices)*1.e-20;



n_dofs_eig (column) = sizeofmatrices ; 

 %apply boundary conditions on all edge nodes
d_indices = [...                       % dofs with constraints
    1:n_nodes_x ...                                     % bottom
 n_nodes_x+1:n_nodes_x:(n_nodes_y-2)*n_nodes_x+1 ...  % left
 2*n_nodes_x:n_nodes_x:(n_nodes_y-1)*n_nodes_x   ...  % right
   (n_nodes_y-1)*n_nodes_x+1:n_nodes...
   
   
   

   ... 1 n_nodes_x n_nodes_x*n_nodes_y n_nodes_x*n_nodes_y-n_nodes_x+1 ...
    ];                  % top


% d_indices=[1 n_nodes_x n_nodes_x+1 n_nodes];

u_indices = setdiff((1:sizeofmatrices),d_indices);
u_indices = setdiff((u_indices),d_indices+n_dofs);
u_indices = setdiff((u_indices),d_indices+n_dofs*2);



J_sub     = jac_mat(u_indices, u_indices);
r_sub     = res_vec(u_indices, u_indices);



V=zeros(sizeofmatrices,4);


 

[V(u_indices,:),D]=eigs(J_sub,r_sub,4,0.1);



% [V,D]=eig(J_sub,r_sub);

D=diag(D);

D(:,2)=1:length(D) ;
D =sortrows(D,1);


% D = D(1:4,:);


D1=(E1*thickness^3)/(12*(1-nu^2));

D(:,1)=sqrt(D(:,1)) * (Lx^2/sqrt(D1/thickness));

 store_eigenvalues(:,column)=D(:,1)
column = column +1; 

if n_elems_x_o  ~= 93

clearvars -except store_eigenvalues column n_dofs_eig reduced  sizeofmatrices V D 
end
end
%  disp_dof = n_dofs;
% 
% n_dofs=(n_dofs_eig(1:4)).^2;
% 
%  exact_num_sol = store_eigenvalues(:,5) ;
% 
% figure(2)
% 
% line1 = abs(store_eigenvalues(1,1:4)-exact_num_sol(1));
% line2 = abs(store_eigenvalues(2,1:4)-exact_num_sol(2));
% line3 = abs(store_eigenvalues(3,1:4)-exact_num_sol(3));
% line4 = abs(store_eigenvalues(4,1:4)-exact_num_sol(4));
% 
% loglog(1./n_dofs,line1,'r-*',1./n_dofs,line2,'g-*',1./n_dofs,line3,'k--*',1./n_dofs,line4,'b-*')
% 
% 
% dof = 1./(n_dofs);
% 
% 
% if reduced == 1
% 
% mat = [dof;line1];
% fileId = fopen('square_mindlin_lagrange_reduced_omega1.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);
% mat = [dof;line2];
% fileId = fopen('square_mindlin_lagrange_reduced_omega2.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);
% mat = [dof;line3];
% fileId = fopen('square_mindlin_lagrange_reduced_omega3.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);
% mat = [dof;line4];
% fileId = fopen('square_mindlin_lagrange_reduced_omega4.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);  
% 
% else 
%     mat = [dof;line1];
% fileId = fopen('square_mindlin_lagrange_full_omega1_thin.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);
% mat = [dof;line2];
% fileId = fopen('square_mindlin_lagrange_full_omega2_thin.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);
% mat = [dof;line3];
% fileId = fopen('square_mindlin_lagrange_full_omega3_thin.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);
% mat = [dof;line4];
% fileId = fopen('square_mindlin_lagrange_full_omega4_thin.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);  
% end
%  
% n_dofs = disp_dof;

% sol=zeros(sizeofmatrices,4);
% 
% sol(:,1)=V(:,D(1,2));
% sol(:,2)=V(:,D(2,2));
% sol(:,3)=V(:,D(3,2));
% sol(:,4)=V(:,D(4,2));
% 
% for i = 1 : 4
% 
%     if sol(2,i) > sol(2+n_nodes_x,i) 
%         
%         sol(:,i)= sol(:,i) * - 1;
%     end
%     storemax_lag(i) = max(sol(1:n_dofs,i));
%     
% end
%-----------------------------------------------------------



% for i = 1:4
%     
%     K(i)=storemax_fit(i) / storemax_lag(i);
%        
%     sol(:,i)=sol(:,i)*K(i);
% end







% 
% % solve the system of equations
% % solve the system
% dX_vec(u_indices)   = J_sub\r_sub;
% % dX_vec   =jac_mat\res_vec;
% % this is the updated solution
% % X_vec    = X_vec + dX_vec;
% X_vec    =  dX_vec;

% X_vec=sol(1:n_nodes,mode);
% % plot 
% 
% j=1;
% o=2;  
% for i_elem_y=1:n_elems_y
%     o=o-1;
% for i_elem_x=1:n_elems_x+1
%     
%     if i_elem_x==n_elems_x/2;
%             for eta=[-1 1]
%                 j=1;
%                for xi=linspace(-1,1,3)
%                    
%                    Ni= 0.25*[(1-xi)*(1-eta); (1+xi)*(1-eta) ; (1+xi)*(1+eta);(1-xi)*(1+eta)];
%                
%                if   xi>=-1  && xi<0 
%                    h=1;
%                    h_=0;
%                end
%                    if xi>=0 && xi<1 
%                        h_=1;
%                      h=0;
%                     end
%                    U_el1=  h*(Ni(1)*X_vec(i_elem_x  +n_nodes_x*(i_elem_y-1))+Ni(2)*X_vec(i_elem_x+1+n_nodes_x*(i_elem_y-1))+ ...
%                        Ni(3)*X_vec(i_elem_x+n_nodes_x+1+n_nodes_x*(i_elem_y-1)  )+Ni(4)*X_vec(i_elem_x+n_nodes_x+n_nodes_x*(i_elem_y-1)))...
%                          +h_*(Ni(1)*X_vec(i_elem_x+2+n_nodes_x*(i_elem_y-1))+Ni(2)*X_vec(i_elem_x+3+n_nodes_x*(i_elem_y-1))+  ...
%                          Ni(3)*X_vec(i_elem_x+n_nodes_x+3+n_nodes_x*(i_elem_y-1))+Ni(4)*X_vec(i_elem_x+n_nodes_x+2+n_nodes_x*(i_elem_y-1)));    
%                    
%                    U_global1(j,o)= U_el1;
%                    
%               
%             j=j+1;
%                end 
%                o=o+1;
%             end
%     
%     end 
% end
% end
% 
% % for i=1:n_elems_y+1
% %     z_vals=x_vec(1:n_nodes/2-1)
% 
% o=0;
% for i_elem_x=1:n_elems_x+1
%     for i_elem_y=1:n_elems_y
%         
%         i_elem                  = (i_elem_y-1)*n_elems_x+i_elem_x;
%         indices_x               = [i_elem_x i_elem_x+1];
%         indices_y               = [i_elem_y i_elem_y+1];
%         indices                 = ...
%             [(indices_y(1)-1)*n_nodes_x+indices_x ...
%             (indices_y(2)-1)*n_nodes_x+fliplr(indices_x)];
% if i_elem_x==n_elems_x/2 || i_elem_x==n_elems_x/2+1  || i_elem_x==n_elems_x/2+2
%    storesol(1,[1:4]+4*o)=indices;
%    storex  (1,[1:2]+2*o)=indices_x;
%   
%    o=o+1;
% end 
%     end
% end 
% o=0;
% p=0;
% a=1:n_nodes;
% for i=1:n_elems_y+1
%     b(1,[1:n_nodes_x/2]+n_nodes_x/2*p)=a([1:n_nodes_x/2]+n_nodes_x/2*(o));
%     o=o+2;
%     p=p+1;
% end
% ele_del=setdiff(b,storesol);
% x_del=setdiff(1:n_nodes_x/2,storex);
% 
% 
% 
% X_ve=X_vec(ele_del);
% z_vals = zeros(length(x_del), n_nodes_y);
% for i_nd=1:n_nodes_y
%     idx            = ((i_nd-1)*length(x_del)+1:i_nd*length(x_del));
%     z_vals(:,i_nd) = X_ve(idx);
% end
% X=x_vals(x_del);
% %surf(x_vals(x_del),y_vals,z_vals')
% %hold on 
% x=linspace(x_vals(n_elems_x/2),x_vals(n_elems_x/2+1),3);
% x=x' ;
% z_vals=z_vals';
% [size_z_vals1 , size_z_vals]=size(z_vals);
% 
% [size_U_global11 , size_U_global1]=size(U_global1');
% 
% z_vals(:,size_z_vals+1:size_z_vals+size_U_global1)=U_global1';
% X       (size_z_vals+1:size_z_vals+size_U_global1)= x;
% surf(X,y_vals,z_vals)
% 
% b=setdiff(1:n_nodes,b);
% ele_del=setdiff(b,storesol);
% x_del=setdiff(n_nodes_x/2:n_nodes_x,storex);
% 
% X_ve=X_vec(ele_del);
% z_val = zeros(length(x_del), n_nodes_y);
% for i_nd=1:n_nodes_y
%     idx            = ((i_nd-1)*length(x_del)+1:i_nd*length(x_del));
%     z_val(:,i_nd) = X_ve(idx);
% end
% 
% X1=x_vals(x_del);
% 
% z_vals(:,size_z_vals+size_U_global1+1:size_z_vals+size_U_global1+size_z_vals)= z_val';
% X       (size_z_vals+size_U_global1+1:size_z_vals+size_U_global1+size_z_vals)= X1;
% figure(5)
% surf(X,y_vals,z_vals)
% 
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('u')
% hold off 
% shading interp
% % title('lagrange method')
% clearvars -except f Lx Ly n_elems_x_o n_elems_y_o d_indices mode storemax_lag store_eigenvalues
%  %% 
% 
% % 
% % n_elems_x = n_elems_x_o+1;
% % n_elems_y = n_elems_y_o;
% % n_nodes_x = n_elems_x+1;
% % n_nodes_y = n_elems_y+1;
% % n_elems   = n_elems_x*n_elems_y;
% % n_nodes   = n_nodes_x*n_nodes_y;
% % n_dofs    = n_nodes;
% % h_x       = Lx / (n_elems_x-1);
% % h_y       = Ly / n_elems_y;
% 
% % mindlin plate theoty 
% % the weak form is given by 
% % int_A(var_kappa' * C_b *kappa)dA  + int_A(var_gamma' * C_s *gamma)dA  +
% % = int_A (var_w * p)dA 
% 
% 
% 
% 
% 
% 
% p = 1;
% E1 = 5.e6 ;
% E2 = 5.e6 ;
% nu=  0.3;
% k = 5/6;
% 
% thickness = Lx/100; 
% 
% 
% C_s1 = ((E1*thickness*k)/(2*(1+nu))) * [1 0 ; 0 1 ];
% 
% C_b1 = ((E1*thickness^3)/(12*(1-nu^2))) * [1 nu 0;nu 1 0;0 0 (1-nu)/2 ];
% 
% C_s2 = ((E2*thickness*k)/(2*(1+nu))) * [1 0 ; 0 1 ];
% 
% C_b2 = ((E2*thickness^3)/(12*(1-nu^2))) * [1 nu 0;nu 1 0;0 0 (1-nu)/2 ];
% 
% n_elems_x = n_elems_x_o+1;
% n_elems_y = n_elems_y_o;
% n_nodes_x = n_elems_x+1;
% n_nodes_y = n_elems_y+1;
% n_elems   = n_elems_x*n_elems_y;
% n_nodes   = n_nodes_x*n_nodes_y;
% n_dofs    = n_nodes;
% h_x       = Lx / n_elems_x;
% h_y       = Ly / n_elems_y;
% 
% % indices for boundary nodes
% bottom_edge_indices   =   (1:n_nodes_x);
% top_edge_indices      =   (n_nodes-n_nodes_x+1:n_nodes);
% left_edge_indices     =   (1:n_nodes_x:n_nodes-n_nodes_x+1);
% right_edge_indices    =   (n_nodes_x:n_nodes_x:n_nodes);
% 
% % nodal locations to plot numerical solution
% x_vals = (0:h_x:Lx)';
% y_vals = (0:h_y:Ly)';
% x_loc  = zeros(n_nodes,2);
% n_elem_nodes = 4;
% for i=1:n_nodes_y
%     x_loc((i-1)*n_nodes_x+1:i*n_nodes_x,1) = x_vals;
%     x_loc((i-1)*n_nodes_x+1:i*n_nodes_x,2) = y_vals(i);
% end
% 
% 
% % four point Gauss-quadrature rule for 2D domain integration
% qp_loc_domain  = ...
%     [-1/sqrt(3)  1/sqrt(3)  -1/sqrt(3)  1/sqrt(3); ...  % xi location
%     -1/sqrt(3)  -1/sqrt(3)  1/sqrt(3)  1/sqrt(3)];     % eta location
% qp_wgt_domain  = [1 1 1 1];
% % two point Gauss-quadrature rule for bottom edge of elem
% qp_loc_bottom  = ...
%     [-1/sqrt(3)  1/sqrt(3); ...  % xi location
%     -1         -1 ];     % eta location
% % two point Gauss-quadrature rule for top edge of elem
% qp_loc_top     = ...
%     [-1/sqrt(3)  1/sqrt(3); ...  % xi location
%     1          1 ];     % eta location
% 
% % two point Gauss-quadrature rule for left edge of elem
% qp_loc_left    = ...
%     [      -1           -1; ...  % xi location
%     -1/sqrt(3)  1/sqrt(3)];      % eta location
% % two point Gauss-quadrature rule for right edge of elem
% qp_loc_right   = ...
%     [       1            1; ...  % xi location
%     -1/sqrt(3)  1/sqrt(3)];      % eta location
% 
% qp_wgt_boundary  = [1 1];
% n_qp_domain  = 4;
% n_qp_boundary= 2;
% 
% 
% jac_mat = zeros(n_dofs*3, n_dofs*3);
% res_vec = zeros(n_dofs*3, n_dofs*3);
% X_vec   = zeros(n_dofs*3, 1);
% dX_vec  = zeros(n_dofs*3, 1);
% indices = zeros(2,1);
% 
% C_s = C_s1;
% 
% C_b = C_b1;
% for i_elem_x=1:n_elems_x
%     
%     if i_elem_x == n_elems_x/2+1
%        C_s = C_s2;
% 
% C_b = C_b2; 
%     end 
%     for i_elem_y=1:n_elems_y
%         i_elem                  = (i_elem_y-1)*n_elems_x+i_elem_x;
%         indices_x               = [i_elem_x i_elem_x+1];
%         indices_y               = [i_elem_y i_elem_y+1];
%         indices                 = ...
%             [(indices_y(1)-1)*n_nodes_x+indices_x ...
%             (indices_y(2)-1)*n_nodes_x+fliplr(indices_x)];
%         
%         
%         w_vec                   = X_vec(indices);   % initial solution estimate
%         tetax                   = X_vec(indices+n_dofs);
%         tetay                   = X_vec(indices+2*n_dofs); 
%         
%         x_vec                   = x_loc(indices,1);
%         y_vec                   = x_loc(indices,2);
%         
%         % zero all matrices before evaluation
%         n_elem_dofs     = max(size(indices));
%         jac_e           = zeros(n_elem_dofs*3, n_elem_dofs*3);
%         res_e           = zeros(n_elem_dofs*3, n_elem_dofs*3);
%         
%         for i_qp = 1:n_qp_domain
%             
%             xi                      = qp_loc_domain(1,i_qp);           % element quadrature point location
%             eta                     = qp_loc_domain(2,i_qp);           % element quadrature point location
%             
%             % solution and derivative quantities at quadrature
%             % point
% [w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
%          , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
% quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);
%             
%             % calculate the residual and Jacobians
%             
%             A = [ zeros(1,4) zeros(1,4) -Bxmat;...
%                   zeros(1,4) Bymat   zeros(1,4);...
%                   zeros(1,4) Bxmat      -Bymat];
%               
%             B = [Bxmat zeros(1,4) -Bmat;...
%                  Bymat  Bmat zeros(1,4)];
%              
%             D = [Bmat zeros(1,4) zeros(1,4)];
%             
%             
%             res_e                   = res_e + qp_wgt_domain(i_qp) * J_det *(...
%                        ...                   % -dphi/dx du/dx
%                       ...                   % -dphi/dy du/dy
%                 + D'  * D);                              % +phi f 
%             
%             jac_e                   = jac_e + qp_wgt_domain(i_qp) * J_det *(...
%                 A' * C_b * A       ...                   
%                +B' * C_s * B );                         
%         end
%         
%         
%         jac_mat(indices, indices)= jac_mat(indices, indices) + jac_e([1:4],[1:4]);
%         
%         
%         jac_mat(indices+n_dofs, indices+n_dofs)= jac_mat(indices+n_dofs, indices+n_dofs) + jac_e([5:8],[5:8]);
%         
%         jac_mat(indices+2*n_dofs, indices+2*n_dofs)= jac_mat(indices+2*n_dofs, indices+2*n_dofs) + jac_e([9:12],[9:12]);
%         
%         
%         jac_mat(indices, indices+n_dofs)= jac_mat(indices, indices+n_dofs) + jac_e([1:4],[5:8]);
%         jac_mat(indices, indices+2*n_dofs)= jac_mat(indices, indices+2*n_dofs) + jac_e([1:4],[9:12]);
%         
%         jac_mat(indices+n_dofs, indices)= jac_mat(indices+n_dofs, indices) + jac_e([5:8],[1:4]);
%         jac_mat(indices+n_dofs, indices+2*n_dofs)= jac_mat(indices+n_dofs, indices+2*n_dofs) + jac_e([5:8],[9:12]);
%         
%         jac_mat(indices+2*n_dofs, indices)= jac_mat(indices+2*n_dofs, indices) + jac_e([9:12],[1:4]);
%         jac_mat(indices+2*n_dofs, indices+n_dofs)= jac_mat(indices+2*n_dofs, indices+n_dofs) + jac_e([9:12],[5:8]);
%          
%         res_vec(indices,indices)         = res_vec(indices,indices) + res_e([1:4],[1:4]);
%         
% %
%     end
% end
% 
% 
% 
% %  apply boundary conditions on all edge nodes
% d_indices = [...                       % dofs with constraints
%     1:n_nodes_x ...                                     % bottom
%  n_nodes_x+1:n_nodes_x:(n_nodes_y-2)*n_nodes_x+1 ...  % left 
%  2*n_nodes_x:n_nodes_x:(n_nodes_y-1)*n_nodes_x   ...  % right 
%     (n_nodes_y-1)*n_nodes_x+1:n_nodes...
%     
%    
%    ];                  % top
% 
% % d_indices=[1 n_nodes_x n_nodes_x+1 n_nodes];
% 
% u_indices = setdiff((1:3*n_nodes),d_indices);
%  u_indices = setdiff(u_indices,n_dofs+d_indices);
%  u_indices = setdiff(u_indices,2*n_dofs+d_indices);
% 
% % u_indices = setdiff((u_indices),[1 n_nodes_x n_nodes-n_elems_x n_nodes]);
% 
% J_sub     = jac_mat(u_indices, u_indices);
% r_sub     = res_vec(u_indices,u_indices);
% 
% 
% V=zeros(3*n_dofs,4);
% 
% 
% [V(u_indices,:),D]=eigs(J_sub,r_sub,4,'smallestabs');
% 
% 
% 
% D=diag(D);
% 
% D(:,2)=1:length(D) ;
% D =sortrows(D,1);
% D1=(E1*thickness^3)/(12*(1-nu^2));
% 
% D(:,1)=sqrt(D(:,1)) * (Lx^2/sqrt(D1/thickness));
% store_eigenvalues_fit(:,1)=D(:,1)
% 
% sol=zeros(3*n_dofs,4);
% 
% sol(:,1)=V(:,D(1,2));
% sol(:,2)=V(:,D(2,2));
% sol(:,3)=V(:,D(3,2));
% sol(:,4)=V(:,D(4,2));
% 
% for i = 1 : 4
% 
%     if sol(2,i) > sol(n_nodes_x+3,i) 
%         
%         sol(:,i)= sol(:,i) * - 1;
%     end
%     storemax_fit(i) = max(sol(1:n_dofs,i));
%     
% end
% %-----------------------------------------------------------
% 
% 
% 
% for i = 1:4
%     
%     K(i)=storemax_lag(i) / storemax_fit(i);
%        
%     sol(:,i)=sol(:,i)*K(i);
% end
% 
% 
% 
% 
% 
% % solve the system of equations
% % solve the system
% dX_vec  = sol(1:n_nodes,mode);
% 
% % this is the updated solution
% % X_vec    = X_vec + dX_vec;
% 
% 
% 
% % plot the solution
% % plot density
% z_vals = zeros(n_nodes_x, n_nodes_y);
% for i_nd=1:n_nodes_y
%     idx            = ((i_nd-1)*n_nodes_x+1:i_nd*n_nodes_x);
%     z_vals(:,i_nd) = dX_vec(idx);
% end
% 
% 
% subplot(1,2,2)
% surf(x_vals, y_vals, z_vals')
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('u')
% title('conventional method')
% 
% hold off
% 

