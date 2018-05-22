function [lag] = interface_terms_calculation(indices,qp_loc,n_elem_dofs,...
    qp_wgt_boundary,n_qp_boundary,n_nodes_x, x_vec, y_vec, w_vec,tetax,tetay,n_dofs, lag,lambda,bending_rigidity,nu)

k1 = 0.5;
k2=0.5;

n1 = [ 1 ; 0];
n2 = [-1 ; 0];

added_nodes = 0 ;
            indices_elm1=indices;
            indices_elm2=[indices(1)+2 indices(2)+2 indices(2)+2+n_nodes_x indices(1)+2+n_nodes_x ];



%----------------------------------

             jac_ee1           = zeros(3*2*n_elem_dofs,3*2*n_elem_dofs);
             jac_ee2           = zeros(3*2*n_elem_dofs,3*2*n_elem_dofs);
             jac_ee3           = zeros(3*2*n_elem_dofs,3*2*n_elem_dofs);
             

 
            for i_qp = 1:n_qp_boundary
                
                xi                      = qp_loc(1,i_qp);           % element quadrature point location
                eta                     = qp_loc(2,i_qp);           % element quadrature point location
 
                % solution and derivative quantities at quadrature
                % point
           
                
                % solution and derivative quantities at quadrature
                % point
                % right edge uses mode = 2
              [w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 2);
                
                % calculate the residual and Jacobians
                  
                  

               
kGh1 = bending_rigidity(2,1);
kGh2 = bending_rigidity(2,2);

D1=bending_rigidity(1,1);
D2=bending_rigidity(1,2);

Q1 = kGh1*[Bxmat zeros(1,4) zeros(1,4) zeros(1,4) Bmat zeros(1,4);Bymat zeros(1,4) -Bmat zeros(1,4) zeros(1,4) zeros(1,4)]; 
Q2 = kGh2*[zeros(1,4)   Bxmat zeros(1,4) zeros(1,4) zeros(1,4) Bmat;zeros(1,4) Bymat zeros(1,4) -Bmat zeros(1,4) zeros(1,4)];

Mxy1 = (D1*(1-nu)/2)* [ zeros(1,4) zeros(1,4) Bxmat zeros(1,4) -Bymat zeros(1,4)];
Mxy2 = (D2*(1-nu)/2)* [ zeros(1,4) zeros(1,4) zeros(1,4) Bxmat zeros(1,4) -Bymat];

w1 = zeros(1,24);
w1(1:4)=Bmat;

w2=zeros(1,24);
w2(5:8)=Bmat;

betay1 = zeros(1,24);
betay1(9:12)= Bmat;

betay2 = zeros(1,24);
betay2(13:16)= Bmat;

betax1 = zeros(1,24);
betax1(17:20)= -Bmat;

betax2 = zeros(1,24);
betax2(21:24)= -Bmat;

jac_ee1 = jac_ee1 +   qp_wgt_boundary(i_qp) * J_det * ( (k1 * Q1' * n1   + k2 * Q2'*n2 )* (w1-w2)  );


jac_ee2 = jac_ee2 +   qp_wgt_boundary(i_qp) * J_det * ( (k1 * Mxy1' * n1(1)   + k2 * Mxy2'*n2(1) )* (betax1-betax2)  );
   
jac_ee3 = jac_ee3 +   qp_wgt_boundary(i_qp) * J_det * ( (k1 * Mxy1' * n1(2)   + k2 * Mxy2'*n2(2) )* (betay1-betay2)  );               
            end
   
            
            
            for i = 1:5
                for j = 1:5
                    
                    jac_ee1(i*4+(1:4),j*4+(1:4)) = jac_ee1(i*4+(1:4),j*4+(1:4));
                    
                    jac_ee2(i*4+(1:4),j*4+(1:4)) = jac_ee2(i*4+(1:4),j*4+(1:4));
                    
                    jac_ee3(i*4+(1:4),j*4+(1:4)) = jac_ee3(i*4+(1:4),j*4+(1:4));
                    
                end
            end
                    
              
            jac_ee=jac_ee1 + jac_ee2 + jac_ee3;

            
            indices = [indices_elm1 indices_elm2 indices_elm1+(n_dofs+added_nodes) indices_elm2+(n_dofs+added_nodes)...
                indices_elm1+2*(n_dofs+added_nodes) indices_elm2+2*(n_dofs+added_nodes)];
            
                        lag(indices,indices)=...
                lag(indices,indices)-jac_ee;
            
            
            
%             jac_e           = zeros(3*2*n_elem_dofs,3*2*n_elem_dofs);
%             jac_e1           = zeros(n_elem_dofs,n_elem_dofs);
%             jac_e2           = zeros(n_elem_dofs,n_elem_dofs);
%             jac_e3           = zeros(n_elem_dofs,n_elem_dofs);
%             jac_e4           = zeros(n_elem_dofs,n_elem_dofs);
% 
%               
%  
%             for i_qp = 1:n_qp_boundary
%                 
%                 xi                      = qp_loc(1,i_qp);           % element quadrature point location
%                 eta                     = qp_loc(2,i_qp);           % element quadrature point location
%  
%                 % solution and derivative quantities at quadrature
%                 % point
%            
%                 
%                 % solution and derivative quantities at quadrature
%                 % point
%                 % right edge uses mode = 2
%                 [u, du_dx, du_dy, Bmat, Bxmat, Bymat, J_det] = ...
%                     quadrature_pt_quantities(xi, eta, x_vec, y_vec, u_vec, mode);
%                 
%                 % calculate the residual and Jacobians
%                 jac_e1           =  jac_e1  + qp_wgt_boundary(i_qp) * J_det *(  [Bxmat;Bymat]'*n1*Bmat    *alph1*k1 ); 
%                 jac_e2           =  jac_e2  + qp_wgt_boundary(i_qp) * J_det *(  [Bxmat;Bymat]'*n2*Bmat    *alph2*k2 ); 
%                  jac_e3          =  jac_e3  + qp_wgt_boundary(i_qp) * J_det *(  [Bxmat;Bymat]'*n1*Bmat    *alph1*k1 );
%                   jac_e4         =  jac_e4  + qp_wgt_boundary(i_qp) * J_det *(  [Bxmat;Bymat]'*n2*Bmat    *alph2*k2 );
%             end
%    
%             ind=1:4;
%             jac_e(ind,ind)=           T'*jac_e1*T;
%             jac_e(ind+4,ind)=         T'*jac_e2*T;
%             jac_e(ind,ind+4)=-        T'*jac_e3*T;
%             jac_e(ind+4,ind+4)=-      T'*jac_e4*T;
%           
%             lag([indices_elm1 indices_elm2],[indices_elm1 indices_elm2])=lag([indices_elm1 indices_elm2],[indices_elm1 indices_elm2])-jac_e;
              
             %% par2
             jac_ee1           = zeros(3*2*n_elem_dofs,3*2*n_elem_dofs);
             jac_ee2           = zeros(3*2*n_elem_dofs,3*2*n_elem_dofs);
             jac_ee3           = zeros(3*2*n_elem_dofs,3*2*n_elem_dofs);
             
 
           
     
                    for i_qp = 1:n_qp_boundary
                    xi                      = qp_loc(1,i_qp);           % element quadrature point location
                    eta                     = qp_loc(2,i_qp);           % element quadrature point location
                    
                    
                   
                    
                    % solution and derivative quantities at quadrature
                    % point
                    % right edge uses mode = 2
              [w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 2);
                    

                    
                    
                    
kGh1 = bending_rigidity(2,1);
kGh2 = bending_rigidity(2,2);

D1=bending_rigidity(1,1);
D2=bending_rigidity(1,2);

Q1 = kGh1*[Bxmat zeros(1,4) zeros(1,4) zeros(1,4) Bmat zeros(1,4);Bymat zeros(1,4) -Bmat zeros(1,4) zeros(1,4) zeros(1,4)]; 
Q2 = kGh2*[zeros(1,4)   Bxmat zeros(1,4) zeros(1,4) zeros(1,4) Bmat;zeros(1,4) Bymat zeros(1,4) -Bmat zeros(1,4) zeros(1,4)];

Mxy1 = (D1*(1-nu)/2)* [ zeros(1,4) zeros(1,4) Bxmat zeros(1,4) -Bymat zeros(1,4)];
Mxy2 = (D2*(1-nu)/2)* [ zeros(1,4) zeros(1,4) zeros(1,4) Bxmat zeros(1,4) -Bymat];

w1 = zeros(1,24);
w1(1:4)=Bmat;

w2=zeros(1,24);
w2(5:8)=Bmat;

betay1 = zeros(1,24);
betay1(9:12)= Bmat;

betay2 = zeros(1,24);
betay2(13:16)= Bmat;

betax1 = zeros(1,24);
betax1(17:20)= -Bmat;

betax2 = zeros(1,24);
betax2(21:24)= -Bmat;

jac_ee1 = jac_ee1 +   qp_wgt_boundary(i_qp) * J_det * (  (w1'-w2') * (k1 * n1' * Q1    +  k2 * n2'*Q2 )  );

jac_ee2 = jac_ee2 +   qp_wgt_boundary(i_qp) * J_det * ( (betax1'-betax2') * (k1 * Mxy1 * n1(1)   + k2 * Mxy2 *n2(1) )  );
   
jac_ee3 = jac_ee3 +   qp_wgt_boundary(i_qp) * J_det * (  (betay1'-betay2')* (k1 * Mxy1  * n1(2)   + k2 * Mxy2 *n2(2) )  ); 
                    
                end
          
       for i = 1:5
                for j = 1:5
                    
                    jac_ee1(i*4+(1:4),j*4+(1:4)) = jac_ee1(i*4+(1:4),j*4+(1:4));
                    
                    jac_ee2(i*4+(1:4),j*4+(1:4)) = jac_ee2(i*4+(1:4),j*4+(1:4));
                    
                    jac_ee3(i*4+(1:4),j*4+(1:4)) = jac_ee3(i*4+(1:4),j*4+(1:4));
                    
                end
       end
            
                
                  jac_ee=jac_ee1 + jac_ee2 + jac_ee3;

            
            indices = [indices_elm1 indices_elm2 indices_elm1+(n_dofs+added_nodes) indices_elm2+(n_dofs+added_nodes)...
                indices_elm1+2*(n_dofs+added_nodes) indices_elm2+2*(n_dofs+added_nodes)];
            
                        lag(indices,indices)=...
                lag(indices,indices)-jac_ee;
            
%             jac_e           = zeros(2*n_elem_dofs,2*n_elem_dofs);
%             jac_e1           = zeros(n_elem_dofs,n_elem_dofs);
%             jac_e2           = zeros(n_elem_dofs,n_elem_dofs);
%             jac_e3           = zeros(n_elem_dofs,n_elem_dofs);
%             jac_e4           = zeros(n_elem_dofs,n_elem_dofs);
%             
% 
%            
%      
%                     for i_qp = 1:n_qp_boundary
%                     xi                      = qp_loc(1,i_qp);           % element quadrature point location
%                     eta                     = qp_loc(2,i_qp);           % element quadrature point location
%                     
%                     
%                    
%                     
%                     % solution and derivative quantities at quadrature
%                     % point
%                     % right edge uses mode = 2
%                     [u, du_dx, du_dy, Bmat, Bxmat, Bymat, J_det] = ...
%                         quadrature_pt_quantities(xi, eta, x_vec, y_vec, u_vec, mode);
%                     
%                     % calculate the residual and Jacobians
%                     jac_e1           =  jac_e1  + qp_wgt_boundary(i_qp) * J_det *(alph1*k1*Bmat'*n1'*[Bxmat;Bymat]);
%                     jac_e2           =  jac_e2  + qp_wgt_boundary(i_qp) * J_det *(alph1*k1*Bmat'*n1'*[Bxmat;Bymat]);
%                     jac_e3           =  jac_e3  + qp_wgt_boundary(i_qp) * J_det *(alph2*k2*Bmat'*n2'*[Bxmat;Bymat]);
%                     jac_e4           =  jac_e4  + qp_wgt_boundary(i_qp) * J_det *(alph2*k2*Bmat'*n2'*[Bxmat;Bymat]);
%                 end
%           
%             ind=1:4;
%             jac_e(ind,ind)=           T'*jac_e1*T;
%             jac_e(ind+4,ind)=-        T'*jac_e2*T;
%             jac_e(ind,ind+4)=+        T'*jac_e3*T;
%             jac_e(ind+4,ind+4)=-      T'*jac_e4*T;
%                 
%              lag([indices_elm1 indices_elm2],[indices_elm1 indices_elm2])=...
%                  lag([indices_elm1 indices_elm2],[indices_elm1 indices_elm2])-jac_e;
   
             %% par3
         
     
%             jac_e           = zeros(2*n_elem_dofs,2*n_elem_dofs);
%             jac_e1           = zeros(n_elem_dofs,n_elem_dofs);
%             jac_e2           = zeros(n_elem_dofs,n_elem_dofs);
%             jac_e3           = zeros(n_elem_dofs,n_elem_dofs);
%             jac_e4           = zeros(n_elem_dofs,n_elem_dofs);
% 
%                
%       
%         
%                     for i_qp = 1:n_qp_boundary
%                     xi                      = qp_loc(1,i_qp);           % element quadrature point location
%                     eta                     = qp_loc(2,i_qp);           % element quadrature point location
%                     
%                     
%             
%                     % solution and derivative quantities at quadrature
%                     % point
%                     % right edge uses mode = 2
%                     [u, du_dx, du_dy, Bmat, Bxmat, Bymat, J_det] = ...
%                         quadrature_pt_quantities(xi, eta, x_vec, y_vec, u_vec, mode);
%                     
%                     % calculate the residual and Jacobians
%                     jac_e1           =  jac_e1  + qp_wgt_boundary(i_qp) * J_det *(lambda*(Bmat'*Bmat));
%                     jac_e2           =  jac_e2  + qp_wgt_boundary(i_qp) * J_det *(lambda*(Bmat'*Bmat));
%                     jac_e3           =  jac_e3  + qp_wgt_boundary(i_qp) * J_det *(lambda*(Bmat'*Bmat));
%                     jac_e4           =  jac_e4  + qp_wgt_boundary(i_qp) * J_det *(lambda*(Bmat'*Bmat));
%                 end
%                 
%             ind=1:4;
%             jac_e(ind,ind)=     T'*jac_e1*T;
%             jac_e(ind+4,ind)=-  T'*jac_e2*T;
%             jac_e(ind,ind+4)=-  T'*jac_e3*T;
%             jac_e(ind+4,ind+4)= T'*jac_e4*T;
%             
%             
% %             mat=zeros(8,8);
% %          
% %             for i=1:8
% %                 mat(i,i)=1.e-3;
% %             end
% %             
%          lag([indices_elm1 indices_elm2],[indices_elm1 indices_elm2])=lag([indices_elm1 indices_elm2],[indices_elm1 indices_elm2])+jac_e;
%          
           
                        jac_ee1           = zeros(2*n_elem_dofs,2*n_elem_dofs);
             jac_ee2           = zeros(2*n_elem_dofs,2*n_elem_dofs);
             jac_ee3           = zeros(2*n_elem_dofs,2*n_elem_dofs);
             
            jac_e1           = zeros(n_elem_dofs,n_elem_dofs);
            jac_e2           = zeros(n_elem_dofs,n_elem_dofs);
            jac_e3           = zeros(n_elem_dofs,n_elem_dofs);
            jac_e4           = zeros(n_elem_dofs,n_elem_dofs);

            jac_e5           = zeros(n_elem_dofs,n_elem_dofs);
            jac_e6           = zeros(n_elem_dofs,n_elem_dofs);
            jac_e7           = zeros(n_elem_dofs,n_elem_dofs);
            jac_e8           = zeros(n_elem_dofs,n_elem_dofs);
      
        
                    for i_qp = 1:n_qp_boundary
                    xi                      = qp_loc(1,i_qp);           % element quadrature point location
                    eta                     = qp_loc(2,i_qp);           % element quadrature point location
                    
                    
            
                    % solution and derivative quantities at quadrature
                    % point
                    % right edge uses mode = 2
              [w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, 2);
                    
                    % calculate the residual and Jacobians
                    jac_e1           =  jac_e1  + qp_wgt_boundary(i_qp) * J_det *(lambda*(Bmat'*Bmat));
                    jac_e2           =  jac_e2  + qp_wgt_boundary(i_qp) * J_det *(lambda*(Bmat'*Bmat));
                    jac_e3           =  jac_e3  + qp_wgt_boundary(i_qp) * J_det *(lambda*(Bmat'*Bmat));
                    jac_e4           =  jac_e4  + qp_wgt_boundary(i_qp) * J_det *(lambda*(Bmat'*Bmat));
                    
                    jac_e5           =  jac_e5  + qp_wgt_boundary(i_qp) * J_det *(lambda*((-Bmat)'*(-Bmat)));
                    jac_e6           =  jac_e6  + qp_wgt_boundary(i_qp) * J_det *(lambda*((-Bmat)'*(-Bmat)));
                    jac_e7           =  jac_e7  + qp_wgt_boundary(i_qp) * J_det *(lambda*((-Bmat)'*(-Bmat)));
                    jac_e8           =  jac_e8  + qp_wgt_boundary(i_qp) * J_det *(lambda*((-Bmat)'*(-Bmat)));
                    
                    
                end
                
            ind=1:4;
            
            jac_ee1(ind,ind)=     jac_e1;
            jac_ee1(ind+4,ind)=-  jac_e2;
            jac_ee1(ind,ind+4)=-  jac_e3;
            jac_ee1(ind+4,ind+4)= jac_e4;

            jac_ee2(ind,ind)=     jac_e1;
            jac_ee2(ind+4,ind)=-  jac_e2;
            jac_ee2(ind,ind+4)=-  jac_e3;
            jac_ee2(ind+4,ind+4)= jac_e4;
            
            jac_ee3(ind,ind)=     jac_e5;
            jac_ee3(ind+4,ind)=-  jac_e6;
            jac_ee3(ind,ind+4)=-  jac_e7;
            jac_ee3(ind+4,ind+4)= jac_e8;
            
%             mat=zeros(8,8);
%          
%             for i=1:8
%                 mat(i,i)=1.e-3;
%             end
%             
         lag([indices_elm1 indices_elm2],[indices_elm1 indices_elm2])=...
             lag([indices_elm1 indices_elm2],[indices_elm1 indices_elm2])+jac_ee1;
         
           
    lag([indices_elm1+(n_dofs+added_nodes) indices_elm2+(n_dofs+added_nodes)],[indices_elm1+(n_dofs+added_nodes) indices_elm2+(n_dofs+added_nodes)])=...
    lag([indices_elm1+(n_dofs+added_nodes) indices_elm2+(n_dofs+added_nodes)],[indices_elm1+(n_dofs+added_nodes) indices_elm2+(n_dofs+added_nodes)])+jac_ee2;            

     lag([indices_elm1+2*(n_dofs+added_nodes) indices_elm2+2*(n_dofs+added_nodes)],[indices_elm1+2*(n_dofs+added_nodes) indices_elm2+2*(n_dofs+added_nodes)])=...
     lag([indices_elm1+2*(n_dofs+added_nodes) indices_elm2+2*(n_dofs+added_nodes)],[indices_elm1+2*(n_dofs+added_nodes) indices_elm2+2*(n_dofs+added_nodes)])+jac_ee3;  
end

