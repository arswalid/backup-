clear 
clc
clf


% p = 1.e3;
E = 5.e6 ;
nu=  0.3;
k = 5/6;
radius=5; % radius
thickness = radius*0.001; 
rho=1;


C_s = ((E*thickness*k)/(2*(1+nu))) * [1 0 ; 0 1 ];

C_b = ((E*thickness^3)/(12*(1-nu^2))) * [1 nu 0;nu 1 0;0 0 (1-nu)/2 ];

mode = 4;


% quadrilateral



a=radius/50; % side a
b=a; % side b
s=[radius-a/2 radius-b/2]; % bottom left start point
% s=[30 -30];
m=[s(1) s(1)+a s(1)+a s(1) s(1);s(2) s(2) s(2)+b s(2)+b s(2)];
plot(m(1,:), m(2,:),'r')
m=m';
hold on
% circle

f = 1.e3;


reduced=1;
column=1;
%square
for n_elems_x = [ 24]% 2^x
n_elems_y = n_elems_x;
n_nodes_x=n_elems_x+1;
n_nodes_y=n_elems_y+1;
n_nodes=n_nodes_x*n_nodes_y;
n_elems_square = n_elems_x*n_elems_y;

h=a/n_elems_x;

x_vec=[s(1) :h :s(1)+a];
y_vec=[s(2) :h :s(2)+a];

x_loc  = zeros(n_nodes,2);
for i=1:n_nodes_y
    x_loc((i-1)*n_nodes_x+1:i*n_nodes_x,1) = x_vec;
    x_loc((i-1)*n_nodes_x+1:i*n_nodes_x,2) = y_vec(i);
end

bottom_edge_indices   =   (1:n_nodes_x);
top_edge_indices      =   (n_nodes-n_nodes_x+1:n_nodes);
left_edge_indices     =   (1:n_nodes_x:n_nodes-n_nodes_x+1);
right_edge_indices    =   (n_nodes_x:n_nodes_x:n_nodes);



iterations=0;
for r=a*(1.2):a:radius % radius
%   C=[0 0]/r;
theta=0:2*pi/360:2*pi; % the angle
C=[radius radius]/r;
m=r*[cos(theta')+C(1) sin(theta')+C(2)]; % the points you asked
l1=plot(m(:,1), m(:,2),'r');
iterations=iterations+1;

end

div=pi/(2*n_elems_x);
numberofdiv = 2*pi /(div);

step=0;
for r=a*(1.2):a:radius % radius
%   C=[0 0]/r;
theta=0+pi/4:div:2*pi+pi/4; % the angle
C=[radius radius]/r;
m=r*[cos(theta')+C(1) sin(theta')+C(2)]; % the points you asked
l1=plot(m(:,1), m(:,2),'r*');
store_ind(1+step:numberofdiv+step,:) = m(1:numberofdiv,:);
step=step+numberofdiv;
end


n_elems_circle =numberofdiv;
n_elems = n_elems_circle*(iterations)+n_elems_square;

store_ind(:,3)=1+n_nodes:length(store_ind)+n_nodes;

bc_indices=store_ind(end-numberofdiv+1:end,3);

element_number=1;

for i_elem_y=1:n_elems_y
    for i_elem_x=1:n_elems_x
      
i_elem                  = (i_elem_y-1)*n_elems_x+i_elem_x;
indices_x               = [i_elem_x i_elem_x+1];
indices_y               = [i_elem_y i_elem_y+1];
indices                 = ...
            [(indices_y(1)-1)*n_nodes_x+indices_x ...
            (indices_y(2)-1)*n_nodes_x+fliplr(indices_x)];
x_vec                   = x_loc(indices,1);
y_vec                   = x_loc(indices,2);

element_number=i_elem;
element_loc_x (1:4,element_number)= [x_vec];
element_loc_y (1:4,element_number)= [y_vec];
element_indice(1:4,element_number) = [indices];


    end
end


% elements_btw=[24 25 26 27; 23 24 27 28;22 23 28 29;21 22 29 30;16 21 30 31;11 16 31 32;6 11 32 33;1 6 33 34;...
%               2 1 34 35;3 2 35 36;4 3 36 37;5 4 37 38;10 5 38 39;15 10 39 40;20 15 40 41;25 20 41 26]';
          
          
elements_btw=[ flip(top_edge_indices(1:end-1)), flip(left_edge_indices(1:end-1)), bottom_edge_indices(2:end), right_edge_indices(2:end);...
                  flip(top_edge_indices(1:end)) flip(left_edge_indices(1:end-1)) bottom_edge_indices(2:end) right_edge_indices(2:end-1);...
                  store_ind(1,3)+[0:numberofdiv-1];...
                 [store_ind(2,3)+[0:numberofdiv-2] store_ind(1,3)]];
          
xloc=[x_loc;store_ind(:,[1 2])];
element_number=element_number+1;
for i=1:numberofdiv
element_indice(1:4,element_number) = elements_btw(:,i);
element_loc_x (1:4,element_number)= [xloc(element_indice(1:4,element_number),1)];
element_loc_y (1:4,element_number)= [xloc(element_indice(1:4,element_number),2)];

element_number=element_number+1;
end



for r = 1:iterations-1
    
    for i = 1:numberofdiv-1
      
element_loc_x (1:4,element_number)= [store_ind(([i+1 i])+(r-1)*numberofdiv,1);store_ind(([i i+1])+(r)*numberofdiv,1)];
element_loc_y (1:4,element_number)= [store_ind(([i+1 i])+(r-1)*numberofdiv,2);store_ind(([i i+1])+(r)*numberofdiv,2)];
element_indice(1:4,element_number) =[store_ind(([i+1 i])+(r-1)*numberofdiv,3);store_ind(([i i+1])+(r)*numberofdiv,3)];
element_number=element_number+1;
    end
element_loc_x (1:4,element_number)= [store_ind(([1 numberofdiv])+(r-1)*numberofdiv,1);store_ind(([numberofdiv 1])+(r)*numberofdiv,1)];
element_loc_y (1:4,element_number)= [store_ind(([1 numberofdiv])+(r-1)*numberofdiv,2);store_ind(([numberofdiv 1])+(r)*numberofdiv,2)];
element_indice(1:4,element_number) =[store_ind(([1 numberofdiv])+(r-1)*numberofdiv,3);store_ind(([numberofdiv 1])+(r)*numberofdiv,3)];
element_number=element_number+1;
end
element_number=element_number-1;


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

n_nodes = max(store_ind(:,3));
n_dofs=n_nodes;

jac_mat = zeros(n_dofs*3, n_dofs*3);
res_vec = zeros(n_dofs*3, n_dofs*3);

jac_mat = sparse(jac_mat);
res_vec =sparse(res_vec);

jac_mat =sparse(jac_mat);
res_vec = sparse (res_vec);

X_vec   = zeros(n_dofs*3, 1);
dX_vec  = zeros(n_dofs*3, 1);
indices = zeros(2,1);
dof(column)=n_dofs*3;
for i_elem =1:element_number
     

        indices                 = element_indice(:,i_elem);
        x_vec                   = element_loc_x(:,i_elem);
        y_vec                   = element_loc_y(:,i_elem);    
        
        
        w_vec                   = X_vec(indices);   % initial solution estimate
        tetax                   = X_vec(indices+n_dofs);
        tetay                   = X_vec(indices+2*n_dofs); 
        

        
        % zero all matrices before evaluation
        n_elem_dofs     = max(size(indices));
        jac_e           = zeros(n_elem_dofs*3, n_elem_dofs*3);
         jac_e1 = zeros(n_elem_dofs*3, n_elem_dofs*3);
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
                     A' * C_b * A  );

        end
        
        if reduced ==1
        for i_qp = 1
            
            xi                      = 0;           % element quadrature point location
            eta                     = 0;           % element quadrature point location
            
            % solution and derivative quantities at quadrature
            % point
            [w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
                , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
                quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 0);
            
            % calculate the residual and Jacobians
            
            
            %
            %             jac_e                   = jac_e + qp_wgt_domain(i_qp) * J_det *(...
            %                 [A;B]' * C * [A;B]       ...
            %                );
            B =  [Bxmat zeros(1,4) -Bmat;...
                 Bymat  Bmat zeros(1,4)];
             
             
                        jac_e1                   =  4 * J_det *(...
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

            B =  [Bxmat zeros(1,4) -Bmat;...
                 Bymat  Bmat zeros(1,4)];
             
             
                        jac_e1                   =  jac_e1 + qp_wgt_domain(i_qp) * J_det *(...
             ...                   
               +B' * C_s * B ); 
                    end
        end
        
     jac_e= jac_e+ jac_e1;

        
        
        
        jac_mat(indices, indices)= jac_mat(indices, indices) + jac_e([1:4],[1:4]);
        
        
        jac_mat(indices+n_dofs, indices+n_dofs)= jac_mat(indices+n_dofs, indices+n_dofs) + jac_e([5:8],[5:8]);
        
        jac_mat(indices+2*n_dofs, indices+2*n_dofs)= jac_mat(indices+2*n_dofs, indices+2*n_dofs) + jac_e([9:12],[9:12]);
        
        
        jac_mat(indices, indices+n_dofs)= jac_mat(indices, indices+n_dofs) + jac_e([1:4],[5:8]);
        jac_mat(indices, indices+2*n_dofs)= jac_mat(indices, indices+2*n_dofs) + jac_e([1:4],[9:12]);
        
        jac_mat(indices+n_dofs, indices)= jac_mat(indices+n_dofs, indices) + jac_e([5:8],[1:4]);
        jac_mat(indices+n_dofs, indices+2*n_dofs)= jac_mat(indices+n_dofs, indices+2*n_dofs) + jac_e([5:8],[9:12]);
        
        jac_mat(indices+2*n_dofs, indices)= jac_mat(indices+2*n_dofs, indices) + jac_e([9:12],[1:4]);
        jac_mat(indices+2*n_dofs, indices+n_dofs)= jac_mat(indices+2*n_dofs, indices+n_dofs) + jac_e([9:12],[5:8]);
         
                res_vec(indices,indices)         = res_vec(indices,indices) + (rho*thickness)* res_e([1:4],1:4);
        res_vec(indices+n_dofs, indices+n_dofs)         = ...
            res_vec(indices+n_dofs, indices+n_dofs) + ((rho*thickness^3)/12)*res_e([1:4],1:4);
        res_vec(indices+2*n_dofs, indices+2*n_dofs)         = ...
            res_vec(indices+2*n_dofs, indices+2*n_dofs) + ((rho*thickness^3)/12)*res_e([1:4],1:4);
%         res_vec(indices,indices)         = res_vec(indices,indices) +  res_e([1:4],1);
end
%     
d_indices=[bc_indices bc_indices+n_dofs bc_indices+2*n_dofs]; %clamped
% d_indices=[bc_indices ]; %simply supported
u_indices=setdiff(1:3*n_dofs,d_indices);

J_sub     = jac_mat(u_indices, u_indices);
r_sub     = res_vec(u_indices, u_indices);



% solve the system of equations
V=zeros(3*n_dofs,10);


[V(u_indices,:),D]=eigs(J_sub,r_sub,10,'smallestabs');



D=diag(D);
D=sqrt(D);

D(:,2)=1:length(D) ;
D =sortrows(D,1);


D1 = (E*thickness^3)/(12*(1-nu^2));

% D(:,1)= D(:,1) * (radius^2) * sqrt(rho/D1);
D(:,1)= D(:,1) * sqrt( (rho*thickness*radius^4)/D1);

store_eigenvalues_fit(:,column)=D(:,1)
column=column+1;

sol=zeros(3*n_dofs,10);

sol(:,1)=V(:,D(1,2));
sol(:,2)=V(:,D(2,2));
sol(:,3)=V(:,D(3,2));
sol(:,4)=V(:,D(4,2));
sol(:,5)=V(:,D(5,2));
sol(:,6)=V(:,D(6,2));
sol(:,7)=V(:,D(7,2));
sol(:,8)=V(:,D(8,2));
sol(:,9)=V(:,D(9,2));
sol(:,10)=V(:,D(10,2));
for i = 1 : 10

    if sol(2,i) > sol(n_nodes_x+3,i) 
        
        sol(:,i)= sol(:,i) * - 1;
    end
    storemax_fit(i) = max(sol(1:n_dofs,i));
    
end


 clearvars -except E nu k radius thickness rho C_s  C_b mode a column store_eigenvalues_fit s reduced  dof
end 

%--------------------convergence study ______
figure(6)
clf
% line1= abs(store_eigenvalues_fit(1,1:4) - store_eigenvalues_fit(1,5))/store_eigenvalues_fit(1,5);
% line2= abs(store_eigenvalues_fit(2,1:4) - store_eigenvalues_fit(2,5))/store_eigenvalues_fit(2,5);
% line3= abs(store_eigenvalues_fit(4,1:4) - store_eigenvalues_fit(4,5))/store_eigenvalues_fit(4,5);
% line4= abs(store_eigenvalues_fit(6,1:4) - store_eigenvalues_fit(6,5))/store_eigenvalues_fit(6,5);

line1= abs(store_eigenvalues_fit(1,1:4) - store_eigenvalues_fit(1,5));
line2= abs(store_eigenvalues_fit(2,1:4) - store_eigenvalues_fit(2,5));
line3= abs(store_eigenvalues_fit(4,1:4) - store_eigenvalues_fit(4,5));
line4= abs(store_eigenvalues_fit(6,1:4) - store_eigenvalues_fit(6,5));

dof=(dof(1:4).^2);
hold on

loglog(1./dof,line1','r-*',1./dof,line2','k-*',1./dof,line3','b-*',1./dof,line4','g-*')

hold on 

legend('1st eigenvalue ','2nd eigenvalue','3rd eigenvalue','4th eigenvalue')
%-----------------------------------------------------------------------

slope = (log(line1(3))-log(line1(2)))/(log(dof(3))-log(dof(2)))

legend('1st','2nd','3rd','4th')
title('full integration')
if reduced==1
%%%%%%%%%%%%%%%%%%%%%%
dof1=1./dof;
mat = [dof1;line1];
fileId = fopen('mindlin_body_reduced_omega1.txt','w');
fprintf (fileId, '%6s %12s\r\n','x','disp');
fprintf (fileId, '%12.8f %12.8f\r\n',mat);
fclose(fileId);
mat = [dof1;line2];
fileId = fopen('mindlin_body_reduced_omega2.txt','w');
fprintf (fileId, '%6s %12s\r\n','x','disp');
fprintf (fileId, '%12.8f %12.8f\r\n',mat);
fclose(fileId);
mat = [dof1;line3];
fileId = fopen('mindlin_body_reduced_omega3.txt','w');
fprintf (fileId, '%6s %12s\r\n','x','disp');
fprintf (fileId, '%12.8f %12.8f\r\n',mat);
fclose(fileId);
mat = [dof1;line4];
fileId = fopen('mindlin_body_reduced_omega4.txt','w');
fprintf (fileId, '%6s %12s\r\n','x','disp');
fprintf (fileId, '%12.8f %12.8f\r\n',mat);
fclose(fileId);
else

    
    %%%%%%%%%%%%%%%%%%%%%%
dof1=1./dof;
mat = [dof1;line1];
fileId = fopen('mindlin_body_full_omega1.txt','w');
fprintf (fileId, '%6s %12s\r\n','x','disp');
fprintf (fileId, '%12.8f %12.8f\r\n',mat);
fclose(fileId);
mat = [dof1;line2];
fileId = fopen('mindlin_body_full_omega2.txt','w');
fprintf (fileId, '%6s %12s\r\n','x','disp');
fprintf (fileId, '%12.8f %12.8f\r\n',mat);
fclose(fileId);
mat = [dof1;line3];
fileId = fopen('mindlin_body_full_omega3.txt','w');
fprintf (fileId, '%6s %12s\r\n','x','disp');
fprintf (fileId, '%12.8f %12.8f\r\n',mat);
fclose(fileId);
mat = [dof1;line4];
fileId = fopen('mindlin_body_full_omega4.txt','w');
fprintf (fileId, '%6s %12s\r\n','x','disp');
fprintf (fileId, '%12.8f %12.8f\r\n',mat);
fclose(fileId);
end