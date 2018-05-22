function [vec_i_el,vec_j_el,z_el] = mat_to_vector(jac_e,ind,ind1,ind2 )

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



vec_i_el = [ind(:,1);ind1(:,1);ind2(:,1);ind(:,1);ind(:,1);ind1(:,1);ind1(:,1);ind2(:,1);ind2(:,1)];
vec_j_el = [ind(:,2);ind1(:,2);ind2(:,2);ind1(:,2);ind2(:,2);ind(:,2);ind2(:,2);ind(:,2);ind1(:,2)];
z_el     = [z1';z5';z9';z4';z7';z2';z8';z3';z6'];


end

