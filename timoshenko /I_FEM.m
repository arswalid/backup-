% code for timoshenko beams using G FEM 
clear 
clc
loop_number = 1;

for n_elems = [ 4 8 16 32 2048 ]-1;
gamma= 1.e13;
Length = 10; 
nu = 0.3;
rho = 1; 
E = 5.e6;
E1 = E;
E2 = E;
kappa = 5/6; 
thickness = Length * 0.1;
I = (thickness^3)/12;
A =  thickness;


lagrange= 1; 
reduced = 1 ; % if 0 ==> full integration 
num_eig = 4 ;
mode = 2;
clamped = 1 ; % if 0 ==> simply supported 




%discretization 
n_nodes = n_elems + 1 ;%first order shape func.
if lagrange == 1
    n_dofs  = 2*n_nodes + 6;
else
    n_dofs  = 2*n_nodes + 4;
end

dofs(loop_number,1)  = n_dofs;


h=Length/n_elems;
x_vec=(0:h:Length)';

            


size_mat1=Length/2;
size_mat2=Length-size_mat1;

n_elem_bef=0;
i=size_mat1;
elemsize = h ;
while i > elemsize
    i=i-elemsize;
    n_elem_bef=n_elem_bef+1;
end
n_elem_dis=n_elem_bef+1;
n_elem_aft=n_elems-n_elem_bef-1;


%boundary
xi=0;
N_boundary = [(1-xi)/2;(1+xi)/2];
vec5 = N_boundary;
dN_boundary = [-.5;.5];
clear xi 

%additional mat properties 
G1 = E1 /(2*(1+nu));
G2 = E2 /(2*(1+nu));
stiff_mat =zeros(n_dofs,n_dofs); 
mass_mat  =zeros(n_dofs,n_dofs); 

    dom =1; left = 0 ; right = 0;
[dN_dN_shear,dN_N,N_dN,N_N_shear,N_N_bend,dN_dN_bend] = elem_mat(reduced,dom,left,right);
if n_elem_bef ~= 0
    E = E1;
    G = G1;
end

mat1 = kappa*G*A * dN_dN_shear *(2/h);
mat5 = E*I * dN_dN_bend * (2/h);
mat2 = kappa*G*A * dN_N ;
mat3 = kappa*G*A *N_dN;
mat4 = kappa*G*A * N_N_shear * (h/2); 
vec1=  rho*A*(h/2) * N_N_bend;
vec2=  rho*I*(h/2) * N_N_bend;

elem_stiff = [mat1,-mat2;-mat3,(mat4+mat5)];
elem_mass  = [vec1,zeros(2,2);zeros(2,2),vec2];

for i =1:n_elem_bef
    
   indices = [i i+1];
     mat_ind = [indices indices+n_nodes+2];
     stiff_mat(mat_ind,mat_ind)=stiff_mat(mat_ind,mat_ind)+elem_stiff;
     mass_mat(mat_ind,mat_ind)= mass_mat(mat_ind,mat_ind) +elem_mass ;
end

if n_elem_aft ~= 0
    E = E2;
    G = G2;
end

mat1 = kappa*G*A * dN_dN_shear *(2/h);
mat5 = E*I * dN_dN_bend * (2/h);
mat2 = kappa*G*A * dN_N ;
mat3 = kappa*G*A *N_dN;
mat4 = kappa*G*A * N_N_shear * (h/2); 
vec1=  rho*A*(h/2) * N_N_bend;
vec2=  rho*I*(h/2) * N_N_bend;

elem_stiff = [mat1,-mat2;-mat3,(mat4+mat5)];
elem_mass  = [vec1,zeros(2,2);zeros(2,2),vec2];
  for i =1:n_elem_aft
     indices = [n_elem_dis+2+i n_elem_dis+3+i];  
     mat_ind = [indices indices+n_nodes+2];
     stiff_mat(mat_ind,mat_ind)=stiff_mat(mat_ind,mat_ind)+elem_stiff;
     mass_mat(mat_ind,mat_ind)= mass_mat(mat_ind,mat_ind) +elem_mass ;
  end
    
% I FEM
i=1;
for element =[n_elem_dis n_elem_dis+1]
   
   if element==n_elem_dis
       indices = [n_elem_dis n_elem_dis+1];
       E=E1;
       G=G1;
      indicesw=indices;
       %------------------------------------------------------
    dom =0; left = 1 ; right = 0;
[dN_dN_shear,dN_N,N_dN,N_N_shear,N_N_bend,dN_dN_bend] = elem_mat(reduced,dom,left,right);
   %------------------------------------------------------
 
   elseif  element==n_elem_dis+1
       E=E2;
       G=G2;
       indices =indices + 2  ; 
       indicesteta = indices;

 dom =0; left = 0 ; right = 1;
[dN_dN_shear,dN_N,N_dN,N_N_shear,N_N_bend,dN_dN_bend] = elem_mat(reduced,dom,left,right);

   end
       
mat1 = kappa*G*A * dN_dN_shear *(2/h);
mat5 = E*I * dN_dN_bend * (2/h);
mat2 = kappa*G*A * dN_N ;
mat3 = kappa*G*A *N_dN;
mat4 = kappa*G*A * N_N_shear * (h/2); 
vec1=  rho*A*(h/2) * N_N_bend;
vec2=  rho*I*(h/2) * N_N_bend;

elem_stiff = [mat1,-mat2;-mat3,(mat4+mat5)];
elem_mass  = [vec1,zeros(2,2);zeros(2,2),vec2];
   
     mat_ind = [indices indices+n_nodes+2];
     stiff_mat(mat_ind,mat_ind)=stiff_mat(mat_ind,mat_ind)+elem_stiff;
     mass_mat(mat_ind,mat_ind)= mass_mat(mat_ind,mat_ind) +elem_mass ;
 
     if lagrange == 1
        stiff_mat(indices,n_dofs-1)=-vec5;
   stiff_mat(n_dofs-1,indices)=-vec5';
      stiff_mat(indices+n_nodes+2,n_dofs)=-vec5;
   stiff_mat(n_dofs,indices+n_nodes+2)=-vec5';
   vec5=-vec5; 
     end
     
   i=i+1;
end


if lagrange == 0
elemsize = h ; 
E_1 = E1;G_1 = G1; k = kappa;
E_2 = E2;G_2 = G2;
[final_mat] = nitsches_timoshenko(N_boundary,dN_boundary,gamma,elemsize,k,G_1,G_2,A,E_1,E_2,I);
   stiff_mat([indices-2 indices indices+n_nodes indices+n_nodes+2],...
       [indices-2 indices indices+n_nodes indices+n_nodes+2])=...
       stiff_mat([indices-2 indices indices+n_nodes indices+n_nodes+2],...
       [indices-2 indices indices+n_nodes indices+n_nodes+2]) + final_mat;
end
   

% bc's
if lagrange == 1 
if clamped 
d_indices= [1 n_nodes+2  n_nodes+3   n_dofs-2  ];%clamped
else
d_indices = [1 n_nodes+2];%simply supported  
end
else 
  if clamped 
d_indices= [1 n_nodes+2  n_nodes+3   n_dofs  ];%clamped
else
d_indices = [1 n_nodes+2];%simply supported  
  end
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
  

if lagrange == 1
solution = V(1:n_dofs-2,mode);
else 
solution = V(1:n_dofs,mode);
end
    
%     plot_IFEM(solution,n_elem_dis,n_nodes,n_elem_bef,n_elem_aft,elemsize)
  
    clearvars -except store_eigenvalues_fitted dofs loop_number 
    
end 
    
    store_eigenvalues_fitted 
    dofs

