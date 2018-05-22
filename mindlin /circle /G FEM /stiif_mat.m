function [jac_e,res_e,n_elem_dofs] = stiif_mat(indices,x_vec, y_vec, w_vec,tetax,tetay,reduced,n_qp_domain,qp_loc_domain,qp_wgt_domain,C_b,C_s )


            % zero all matrices before evaluation
            n_elem_dofs     = max(size(indices));
            jac_e           = zeros(n_elem_dofs*3,n_elem_dofs*3);
       
            res_e           = zeros(n_elem_dofs*3,n_elem_dofs*3);
            
            
            for i_qp = 1:n_qp_domain
                
                xi                      = qp_loc_domain(1,i_qp);           % element quadrature point location
                eta                     = qp_loc_domain(2,i_qp);           % element quadrature point location
                
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
               B' * C_s * B );       
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
               B' * C_s * B );       
                                        end
        end
end

