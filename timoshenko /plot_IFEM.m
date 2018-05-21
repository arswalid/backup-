function [] = plot_IFEM(solution,n_elem_dis,n_nodes,n_elem_bef,n_elem_aft,elemsize)
clf



    if solution(2) > solution(3) 
        solution= solution * - 1;
    end
    
%-----------------------------------------------------------


    
    K=1 / max(solution);
    solution=solution*K;

i=1;
n=3;


u = solution;

for xi = linspace(-1,1,n)
     
   
    if   xi>=-1  && xi<0 
                   h=1;
                   h_=0;
               
    elseif xi>=0 && xi<1 
                     h_=1;
                     h=0;
                    end
    
    N = [(1-xi)/2;(1+xi)/2];
    
    U_global  (i)= h*(N(1)*solution(n_elem_dis)+N(2)*solution(n_elem_dis+1)) +...
        h_*(N(1)*solution(n_elem_dis+2)+N(2)*solution(n_elem_dis+3));
    
    rot_glob  (i)= h*(N(1)*solution(n_elem_dis+n_nodes+2)+N(2)*solution(n_elem_dis+n_nodes+2+1))...
        + h_*(N(1)*solution(n_elem_dis+n_nodes+2+2)+N(2)*solution(n_elem_dis+n_nodes+2+3));
    i=i+1;
    
end
            h=elemsize;
            x=linspace(h*n_elem_bef,h*(n_elem_bef+1),n);
     k=1;
     
while n_elem_bef>0
U(k)=u(k);
R(k)=u(k+n_nodes+2);
X(k)=k*h-h;
k=k+1;
n_elem_bef=n_elem_bef-1;
end

U(k:n+k-1)=U_global;
R(k:n+k-1)=rot_glob;
i=k;
X(k:n+k-1)=x;


while n_elem_aft>0
U(n+k)=u(k+4);
R(n+k)=u(k+n_nodes+2+4);
X(n+k)=k*h+h;

k=k+1;
n_elem_aft=n_elem_aft-1;
end  

subplot(2,1,1)
% plot([X(1:i) X(n-1+i:end)],[U(1:i) U(n-1+i:end)],'ro')
            hold on         %dof plot
            subplot(2,1,2)
% plot([X(1:i) X(n-1+i:end)],[R(1:i) R(n-1+i:end)],'ro')
            hold on         %dof plot
            subplot(2,1,1)
            plot(X,U,'k')
%             xticks([0:h:L])
            hold on 
            subplot(2,1,2)
            plot(X,R,'k')
            hold on 
end

