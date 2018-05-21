function [dN_dN_shear,dN_N,N_dN,N_N_shear,N_N_bend,dN_dN_bend] = elem_mat(reduced,dom,left,right)



if reduced ==1 
   dN_dN_shear=zeros(2,2);
    dN_N=zeros(2,2);
    N_dN=zeros(2,2);
    N_N_shear=zeros(2,2);
   
    if dom == 1 
    xi =  0;
    w  =  2;
    elseif right == 1 
    xi = 0.5;
    w  =  1 ;
    elseif left == 1
    xi = -0.5;
    w  = 1;
    end
    
    for i = 1
        
    dN_dN_shear  = dN_dN_shear+ w(i) * [1/4 -1/4;-1/4 1/4];
    
    dN_N   =dN_N+ w(i) * [(conj(xi(i))/4 -1/4), (-conj(xi(i))/4 -1/4);...
        (1/4 -conj(xi(i))/4), ( conj(xi(i))/4 +1/4)];
    
    N_dN  =N_dN+ w(i) * [( xi(i)/4 - 1/4), (1/4 -xi(i)/4);...
        (-xi(i)/4 - 1/4), (xi(i)/4 +1/4)];
    
    N_N_shear = N_N_shear  + w(i) * [((xi(i)/2 -1/2)*(conj(xi(i))/2 -1/2)) , (-(xi(i)/2 - 1/2)*(conj(xi(i))/2 + 1/2));...
        (-(xi(i)/2 + 1/2)*(conj(xi(i))/2 -1/2)), ((xi(i)/2 +1/2)*(conj(xi(i))/2 + 1/2))];
    end

    
elseif reduced == 0 
    dN_dN_shear=zeros(2,2);
    dN_N=zeros(2,2);
    N_dN=zeros(2,2);
    N_N_shear=zeros(2,2);

    if dom == 1 
    xi = [-1/sqrt(3)  1/sqrt(3)];
    w  = [ 1 1];
    elseif right == 1
    xi = [0.2113 0.7887];
    w  = [0.5 0.5 ];
    elseif left ==1
    xi = [-0.7887 -0.2113];
    w  = [0.5 0.5 ];
    end
    

    
    for i = 1:2
        
    dN_dN_shear  = dN_dN_shear+ w(i) * [1/4 -1/4;-1/4 1/4];
    
    dN_N   =dN_N+ w(i) * [(conj(xi(i))/4 -1/4), (-conj(xi(i))/4 -1/4);...
        (1/4 -conj(xi(i))/4), ( conj(xi(i))/4 +1/4)];
    
    N_dN  =N_dN+ w(i) * [( xi(i)/4 - 1/4), (1/4 -xi(i)/4);...
        (-xi(i)/4 - 1/4), (xi(i)/4 +1/4)];
    
    N_N_shear = N_N_shear  + w(i) * [((xi(i)/2 -1/2)*(conj(xi(i))/2 -1/2)) , (-(xi(i)/2 - 1/2)*(conj(xi(i))/2 + 1/2));...
        (-(xi(i)/2 + 1/2)*(conj(xi(i))/2 -1/2)), ((xi(i)/2 +1/2)*(conj(xi(i))/2 + 1/2))];
    end
    
end
    
    % 2 point quadrature weight = 1 loc = +- 1/sqrt(3);

    if dom == 1 
    xi = [-1/sqrt(3)  1/sqrt(3)];
    w  = [ 1 1];
    elseif right ==1
    xi = [0.2113 0.7887];
    w  = [0.5 0.5 ];
    elseif left ==1
    xi = [-0.7887 -0.2113];
    w  = [0.5 0.5 ];
    end
    
    N_N_bend=zeros(2,2);
    dN_dN_bend =zeros(2,2);
    
    for i=1:2
        
        N_N_bend = N_N_bend + w(i)*([(1-xi(i))/2;(1+xi(i))/2]*[(1-xi(i))/2;(1+xi(i))/2]');
        
        dN_dN_bend  = dN_dN_bend +  w(i) * [1/4 -1/4 ;-1/4 1/4];
    end
    
    
    
end