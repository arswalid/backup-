% code for timoshenko beams using G FEM 
clear 
clc
loop_number = 1;
for n_elems = [ 16 ]
Length = 10; 
nu = 0.3;
rho = 1; 
E = 5.e6;
E1 = E;
E2 = E;
kappa = 5/6; 
thickness = Length * 0.01;
I = (thickness^3)/12;
A =  thickness;



reduced = 0 ; % if 0 ==> full integration 
num_eig = 4 ;
mode = 1;
clamped = 1 ; % if 0 ==> simply supported 




%discretization 
n_nodes = n_elems + 1 ;%first order shape func.
h = Length/n_elems;
n_dofs = n_nodes * 2; % disp and rotation 
x_vec = (0:h:Length)';

dofs(loop_number,1) = n_dofs ; 

stiff_mat =zeros(n_nodes*2,n_nodes*2 );
mass_mat =zeros(n_nodes*2,n_nodes*2);

%additional mat properties 
G1 = E1 /(2*(1+nu));
G2 = E2 /(2*(1+nu));

%material 1 
dom = 1; left =0; right =0;
[dN_dN_shear,dN_N,N_dN,N_N_shear,N_N_bend,dN_dN_bend] = elem_mat(reduced,dom,left,right);



for i =1:n_elems
    
 if i <=   n_elems/2 
E = E1;
G = G1;
 else 
E = E2;
G = G2;
 end

mat1 = kappa*G*A * dN_dN_shear *(2/h);
mat2 = kappa*G*A * dN_N ;
mat3 = kappa*G*A * N_dN;
mat4 = kappa*G*A * N_N_shear * (h/2); 
mat5 = E*I * dN_dN_bend * (2/h);
vec1=  rho*A* (h/2) * N_N_bend;
vec2=  rho*I*(h/2) * N_N_bend;

elem_stiff = [mat1,-mat2;-mat3,(mat4+mat5)];
elem_mass  = [vec1,zeros(2,2);zeros(2,2),vec2];
    
    
     indices = [i i+1]; 
     mat_ind = [indices indices+n_nodes];
     stiff_mat(mat_ind,mat_ind)=stiff_mat(mat_ind,mat_ind)+elem_stiff;
     mass_mat(mat_ind,mat_ind)= mass_mat(mat_ind,mat_ind) +elem_mass ;
end

% bc's
if clamped 
d_indices = [1 n_nodes n_nodes+1 n_dofs];%clamped
else
d_indices = [1 n_nodes];%simply supported  
end

u_indices = setdiff(1:n_dofs,d_indices);


stiff_mat_sub = stiff_mat(u_indices,u_indices);
mass_mat_sub  = mass_mat(u_indices,u_indices);

 V=zeros(n_dofs,num_eig);
[V(u_indices,:),D]=eigs(stiff_mat_sub,mass_mat_sub,num_eig,'smallestabs');

D=diag(D);
% store_eigenvalues_fitted(:,loop_number)=D
store_eigenvalues_fitted(:,loop_number)=sqrt(sqrt(D) * Length^2 *sqrt(rho *A/(E*I)) );
loop_number=loop_number+1;

% plot __________________________
clf
solution = V(1:n_nodes,mode);
if solution(2) > solution(3)
    solution = solution * -1 ; 
end
%scaling 
scale_fac = 1/max(solution) ; 
solution = scale_fac * solution ; 
plot(x_vec,solution,'--')
hold on
%________________________________

clearvars -except store_eigenvalues_fitted dofs loop_number
end
store_eigenvalues_fitted
dofs

