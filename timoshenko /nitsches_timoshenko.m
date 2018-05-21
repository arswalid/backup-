function [final_mat] = nitsches_timoshenko(N_boundary,dN_boundary,gamma,elemsize,k,G_1,G_2,A,E_1,E_2,I)



% mat1 = N_boundary * dN_boun dary' * 0.5 * (2/elemsize);
% 
% term3([1 2 ],[1 2]) =  mat1; 
% term3([3 4 ],[1 2]) = -mat1; 
% term3([1 2 ],[3 4]) = -mat1; 
% term3([3 4],[3 4])  =  mat1; 
% 
% term3([5:8],[5:8])=term3([1:4],[1:4]);
% 
% 
% mat2 = dN_boundary * N_boundary' * 0.5 * (2/elemsize);
% 
% term2([1 2 ],[1 2]) =  mat2; 
% term2([3 4 ],[1 2]) = -mat2; 
% term2([1 2 ],[3 4]) = -mat2; 
% term2([3 4],[3 4])  =  mat2; 
% 
% term2([5:8],[5:8])=term2([1:4],[1:4]);

term1=zeros(8,8);

mat1 =  N_boundary *  N_boundary' * 0.5 *G_1*A*k* (2/elemsize);
mat2 = -N_boundary * dN_boundary' * 0.5 *G_1*A*k* (2/elemsize);
mat3 = -N_boundary *  N_boundary' * 0.5 *G_2*A*k* (2/elemsize);
mat4 =  N_boundary * dN_boundary' * 0.5 *G_2*A*k* (2/elemsize);

term1(1:2,5:6) = term1(1:2,5:6)+mat1;
term1(1:2,1:2) = term1(1:2,1:2)+mat2;
term1(1:2,7:8) = term1(1:2,7:8)+mat3;
term1(1:2,3:4) = term1(1:2,3:4)+mat4;

term1(3:4,5:6) = term1(3:4,5:6)-mat1;
term1(3:4,1:2) = term1(3:4,1:2)-mat2;
term1(3:4,7:8) = term1(3:4,7:8)-mat3;
term1(3:4,3:4) = term1(3:4,3:4)-mat4;

term2=zeros(8,8);

mat1 =  N_boundary *  N_boundary' * 0.5 *G_1*A*k* (2/elemsize);
mat2 = -dN_boundary *  N_boundary' * 0.5 *G_1*A*k* (2/elemsize);
mat3 = -N_boundary *  N_boundary' * 0.5 *G_2*A*k* (2/elemsize);
mat4 =  dN_boundary *  N_boundary' * 0.5 *G_2*A*k* (2/elemsize);

term2(5:6,1:2) = term2(5:6,1:2)+mat1;
term2(1:2,1:2) = term2(1:2,1:2)+mat2;
term2(7:8,1:2) = term2(7:8,1:2)+mat3;
term2(3:4,1:2) = term2(3:4,1:2)+mat4;

term2(5:6,3:4) = term2(5:6,3:4)-mat1;
term2(1:2,3:4) = term2(1:2,3:4)-mat2;
term2(7:8,3:4) = term2(7:8,3:4)-mat3;
term2(3:4,3:4) = term2(3:4,3:4)-mat4;


mat1 = N_boundary * dN_boundary' * 0.5 * E_1*I*(2/elemsize);
mat2 = N_boundary * dN_boundary' * 0.5 * E_2*I*(2/elemsize);

term3([1 2 ],[1 2]) =  mat1; 
term3([3 4 ],[1 2]) = -mat1; 
term3([1 2 ],[3 4]) = -mat2; 
term3([3 4] ,[3 4])  =  mat2; 

term3([5:8],[5:8])=term3([1:4],[1:4]);
term3([1:4],[1:4])=zeros(4,4);


mat2 = dN_boundary * N_boundary' * 0.5 * E_1*I*(2/elemsize);
mat1 = dN_boundary * N_boundary' * 0.5 * E_2*I*(2/elemsize);

term4([1 2 ],[1 2]) =  mat2; 
term4([3 4 ],[1 2]) = -mat1; 
term4([1 2 ],[3 4]) = -mat2; 
term4([3 4],[3 4])  =  mat1; 

term4([5:8],[5:8])=term4([1:4],[1:4]);
term4([1:4],[1:4])=zeros(4,4);





   augmented =zeros(8,8);
   
   term5 = N_boundary * N_boundary';
   
   augmented([1 2 ],[1 2]) = term5;
   augmented([3 4 ],[1 2]) = -term5;
   augmented([1 2 ],[3 4]) = -term5;
   augmented([3 4 ],[3 4]) = term5;
   
   augmented([5 6 ],[5 6]) = term5;
   augmented([7 8 ],[5 6]) = -term5;
   augmented([5 6 ],[7 8]) = -term5;
   augmented([7 8 ],[7 8]) = term5;
   
   augmented=augmented*gamma;


final_mat = augmented - term2 -term3 - term1-term4;
end

