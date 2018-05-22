% mindlin plate theoty 
% the weak form is given by 
% int_A(var_kappa' * C_b *kappa)dA  + int_A(var_gamma' * C_s *gamma)dA  +
% = int_A (var_w * p)dA 

%test

clear 
clc
clf 


E = 5.e6 ;
nu=  0.3;
k = 5/6;
Lx= 10;
Ly= 10;
thickness = Lx*0.01; 
rho  = 1;
mode = 8;

C_s = ((E*thickness*k)/(2*(1+nu))) * [1 0 ; 0 1 ];

C_b = ((E*thickness^3)/(12*(1-nu^2))) * [1 nu 0;nu 1 0;0 0 (1-nu)/2 ];

step=1;
% [  35 45 55 93 149 159]+1;
dofs=[35];
numb_var = length(dofs);
store_eigenvalues_fit=zeros(10,numb_var);
n_dofs_eig=zeros(numb_var);
for n_elems_x = dofs
    
    
    
n_elems_y = n_elems_x;

reduced = 1;


n_nodes_x = n_elems_x+1;
n_nodes_y = n_elems_y+1;
n_elems   = n_elems_x*n_elems_y;
n_nodes   = n_nodes_x*n_nodes_y;
n_dofs    = n_nodes;
h_x       = Lx / n_elems_x;
h_y       = Ly / n_elems_y;



n_dofs_eig(step)=n_dofs*3;


% indices for boundary nodes
bottom_edge_indices   =   (1:n_nodes_x);
top_edge_indices      =   (n_nodes-n_nodes_x+1:n_nodes);
left_edge_indices     =   (1:n_nodes_x:n_nodes-n_nodes_x+1);
right_edge_indices    =   (n_nodes_x:n_nodes_x:n_nodes);

% nodal locations to plot numerical solution
x_vals = (0:h_x:Lx)';
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
% two point Gauss-quadrature rule for top edge of elem
qp_loc_top     = ...
    [-1/sqrt(3)  1/sqrt(3); ...  % xi location
    1          1 ];     % eta location
% two point Gauss-quadrature rule for left edge of elem
qp_loc_left    = ...
    [      -1           -1; ...  % xi location
    -1/sqrt(3)  1/sqrt(3)];      % eta location
% two point Gauss-quadrature rule for right edge of elem
qp_loc_right   = ...
    [       1            1; ...  % xi location
    -1/sqrt(3)  1/sqrt(3)];      % eta location

qp_wgt_boundary  = [1 1];
n_qp_domain  = 4;
n_qp_boundary= 2;


jac_mat = sparse(n_dofs*3, n_dofs*3);
res_vec = sparse(n_dofs*3, n_dofs*3);
vec_i_stiff=0;
vec_j_stiff=0;
z_stiff=0;
vec_i_mass=0;
vec_j_mass=0;
z_mass=0;

elem1=1;

X_vec   = zeros(n_dofs*3, 1);
dX_vec  = zeros(n_dofs*3, 1);
indices = zeros(2,1);
i=1;
for i_elem_x=1:n_elems_x
    for i_elem_y=1:n_elems_y
        i_elem                  = (i_elem_y-1)*n_elems_x+i_elem_x;
        indices_x               = [i_elem_x i_elem_x+1];
        indices_y               = [i_elem_y i_elem_y+1];
        indices                 = ...
            [(indices_y(1)-1)*n_nodes_x+indices_x ...
            (indices_y(2)-1)*n_nodes_x+fliplr(indices_x)];
        
        
        w_vec                   = X_vec(indices);   % initial solution estimate
        tetax                   = X_vec(indices+n_dofs);
        tetay                   = X_vec(indices+2*n_dofs); 
        
        x_vec                   = x_loc(indices,1);
        y_vec                   = x_loc(indices,2);
        
        
        if i == 1
            
            % zero all matrices before evaluation
            n_elem_dofs     = max(size(indices));
            jac_e           = zeros(n_elem_dofs*3, n_elem_dofs*3);
            res_e           = zeros(n_elem_dofs*3, n_elem_dofs*3);
            
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
                    
                    % calculate the residual and Jacobians
                    
                    B =  [Bxmat zeros(1,4) -Bmat;...
                        Bymat  Bmat zeros(1,4)];
                    
                    
                    jac_e                   = jac_e + 4 * J_det *(...
                        ...
                        +B' * C_s * B );
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
                    
                    % calculate the residual and Jacobians
                    
                    
                    
                    B = [Bxmat zeros(1,4) -Bmat;...
                        Bymat  Bmat zeros(1,4)];
                    
                    
                    jac_e                   = jac_e + qp_wgt_domain(i_qp) * J_det *(...
                        ...
                        +B' * C_s * B );
                    
                end
                
            end
            
        end
        
        i=i+1;
        
        
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


jac_mat = sparse(vec_i_stiff,vec_j_stiff,z_stiff,3*n_dofs,3*n_dofs);
res_vec = sparse(vec_i_mass,vec_j_mass,z_mass,3*n_dofs,3*n_dofs);

%  apply boundary conditions on all edge nodes
d_indices = [...                       % dofs with constraints
    1:n_nodes_x ...                                     % bottom
  ...n_nodes_x+1:n_nodes_x:(n_nodes_y-2)*n_nodes_x+1 ...  % left 
   ... 2*n_nodes_x:n_nodes_x:(n_nodes_y-1)*n_nodes_x   ...  % right 
   (n_nodes_y-1)*n_nodes_x+1:n_nodes
   ];                  % top

% d_indices1=[1 n_nodes_x n_nodes_x+1 n_nodes];

u_indices = setdiff((1:3*n_nodes),d_indices);
% u_indices = setdiff(u_indices,d_indices1);
u_indices = setdiff(u_indices,n_dofs+d_indices);
% u_indices = setdiff(u_indices,d_indices1+n_dofs);
% u_indices = setdiff(u_indices,d_indices1+2*n_dofs);
u_indices = setdiff(u_indices,2*n_dofs+d_indices);

% u_indices = setdiff((u_indices),[1 n_nodes_x n_nodes-n_elems_x n_nodes]);

J_sub     = jac_mat(u_indices, u_indices);
r_sub     = res_vec(u_indices, u_indices);



% solve the system of equations
V=zeros(3*n_dofs,10);


[V(u_indices,:),D]=eigs(J_sub,r_sub,10,'smallestabs');



D=diag(D)

D1=(E*thickness^3)/(12*(1-nu^2));

D(:,1)=sqrt(D) * (Lx^2/sqrt(D1/thickness));


store_eigenvalues_fit(:,step)=D(:,1);
step=step+1;

sol=zeros(3*n_dofs,10);

sol(:,1)=V(:,1);
sol(:,2)=V(:,2);
sol(:,3)=V(:,3);
sol(:,4)=V(:,4);
sol(:,5)=V(:,5);
sol(:,6)=V(:,6);
sol(:,7)=V(:,7);
sol(:,8)=V(:,8);
sol(:,9)=V(:,9);
sol(:,10)=V(:,10);

for i = 1 : 10

    if sol(2,i) > sol(n_nodes_x+3,i) 
        
        sol(:,i)= sol(:,i) * - 1;
    end
    storemax_fit(i) = max(sol(1:n_dofs,i));
    
end
%-----------------------------------------------------------



% solve the system of equations
% solve the system
dX_vec  = sol(1:n_nodes,mode);
X_vec=dX_vec ;



% plot the solution
% plot density
z_vals = zeros(n_nodes_x, n_nodes_y);
for i_nd=1:n_nodes_y
    idx            = ((i_nd-1)*n_nodes_x+1:i_nd*n_nodes_x);
    z_vals(:,i_nd) = X_vec(idx);
end
figure(1)
surf(x_vals, y_vals, z_vals')
xlabel('x (m)')
ylabel('y (m)')
zlabel('u')

clearvars -except store_eigenvalues_fit step E nu k Lx Ly thickness mode C_s C_b  n_dofs_eig rho thickness dofs 
end

store_eigenvalues_fit

% plot(n_dofs_eig,store_eigenvalues_fit(1,1:5))

% n_dofs=(n_dofs_eig(1:4)).^2;
% 
%  exact_num_sol = store_eigenvalues_fit(:,5);
% 
% figure(2)
% clf
% 
% line1 = abs(store_eigenvalues_fit(1,1:4)-exact_num_sol(1));
% line2 = abs(store_eigenvalues_fit(2,1:4)-exact_num_sol(2));
% line3 = abs(store_eigenvalues_fit(3,1:4)-exact_num_sol(3));
% line4 = abs(store_eigenvalues_fit(4,1:4)-exact_num_sol(4));
% 
% % line1 = abs(store_eigenvalues_fit(1,:)-exact_num_sol(1));
% % line2 = abs(store_eigenvalues_fit(2,:)-exact_num_sol(2));
% % line3 = abs(store_eigenvalues_fit(3,:)-exact_num_sol(3));
% % line4 = abs(store_eigenvalues_fit(4,:)-exact_num_sol(4));
% 
% loglog(1./n_dofs,line1,'r-*',1./n_dofs,line2,'g-*',1./n_dofs,line3,'k--*',1./n_dofs,line4,'b-*')
% 
% slope1= (log(line1(4))-log(line1(3)))/(log(n_dofs(4))-log(n_dofs(3)));
% slope2= (log(line1(3))-log(line1(2)))/(log(n_dofs(3))-log(n_dofs(2)));
% slope3= (log(line1(2))-log(line1(1)))/(log(n_dofs(2))-log(n_dofs(1)));
% 
% slope11=(slope1+slope2+slope3)/3;
% 
% slope1= (log(line2(4))-log(line2(3)))/(log(n_dofs(4))-log(n_dofs(3)));
% slope2= (log(line2(3))-log(line2(2)))/(log(n_dofs(3))-log(n_dofs(2)));
% slope3= (log(line2(2))-log(line2(1)))/(log(n_dofs(2))-log(n_dofs(1)));
% 
% slope22=(slope1+slope2+slope3)/3;
% 
% slope1= (log(line3(4))-log(line3(3)))/(log(n_dofs(4))-log(n_dofs(3)));
% slope2= (log(line3(3))-log(line3(2)))/(log(n_dofs(3))-log(n_dofs(2)));
% slope3= (log(line3(2))-log(line3(1)))/(log(n_dofs(2))-log(n_dofs(1)));
% 
% slope33=(slope1+slope2+slope3)/3;
% 
% slope1= (log(line4(4))-log(line4(3)))/(log(n_dofs(4))-log(n_dofs(3)));
% slope2= (log(line4(3))-log(line4(2)))/(log(n_dofs(3))-log(n_dofs(2)));
% slope3= (log(line4(2))-log(line4(1)))/(log(n_dofs(2))-log(n_dofs(1)));
% 
% slope44=(slope1+slope2+slope3)/3;
% 
% slope=(slope11+slope22+slope33+slope44)/4
% 
% 
% double(line4(1))*100
% double(line1(4))*100
% 
% dof = 1./n_dofs;
% mat = [dof;line1];
% fileId = fopen('square_mindlin_body_full_omega1.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);
% mat = [dof;line2];
% fileId = fopen('square_mindlin_body_full_omega2.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);
% mat = [dof;line3];
% fileId = fopen('square_mindlin_body_full_omega3.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);
% mat = [dof;line4];
% fileId = fopen('square_mindlin_body_full_omega4.txt','w');
% fprintf (fileId, '%6s %12s\r\n','x','disp');
% fprintf (fileId, '%12.8f %12.8f\r\n',mat);
% fclose(fileId);  
% 
% legend('1st eigenvalue','2nd eigenvalue','3rd eigenvalue','4th eigenvalue')
