function [w , dw_dx , dw_dy , betax , dbetax_dx , dbetax_dy , betay ...
         , dbetay_dx , dbetay_dy , Bmat, Bxmat, Bymat, J_det] = ...
quadrature_pt_quantities(xi, eta, x_vec, y_vec, w_vec,tetax,tetay, mode)


Ni= 0.25*[(1-xi)*(1-eta); (1+xi)*(1-eta) ; (1+xi)*(1+eta);(1-xi)*(1+eta)];
dNi_dxi=0.25*[-(1-eta); (1-eta) ; (1+eta);-(1+eta)];
dNi_deta=0.25*[-(1-xi); -(1+xi) ; (1+xi);(1-xi)];
dx_dxi=dNi_dxi'*x_vec;
dx_deta=dNi_deta'*x_vec;
dy_dxi=dNi_dxi'*y_vec;
dy_deta=dNi_deta'*y_vec;
J=[dx_dxi  dx_deta;dy_dxi dy_deta];
J_inv=inv(J');
dNi_dx=J_inv(1,:)*[dNi_dxi';dNi_deta'];
dNi_dy=J_inv(2,:)*[dNi_dxi';dNi_deta'];
Bmat=Ni';
Bxmat=dNi_dx;
Bymat=dNi_dy;
w=Bmat*w_vec;
betax =-Bmat*tetay;
betay = Bmat*tetax;


dw_dx=Bxmat*w_vec;
dw_dy=Bymat*w_vec;

dbetax_dx=-Bxmat*tetay;
dbetax_dy=-Bymat*tetay;

dbetay_dx=Bxmat*tetax;
dbetay_dy=Bymat*tetax;

if mode==0
    J_det=det(J);
elseif mode==1 || mode==3
    vec1=J';
    vec2= vec1(1,:)';
    J_det= sqrt(vec2'*vec2);
elseif mode==2 || mode==4
    vec1=J';
    vec2= vec1(2,:)';
    J_det= sqrt(vec2'*vec2);
end
