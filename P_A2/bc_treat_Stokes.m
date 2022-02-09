function [a_S, f_S] =bc_treat_Stokes(nx,ny,node_xy,a_S_linear,f,t)



node_indx_u=[1:(2*nx-1)*(2*ny-1)];
BC_node1_u=union([1:(2*nx-1)], [(2*nx-1)*(2*ny-1)-(2*nx-1)+1:(2*nx-1)*(2*ny-1)]);   % lower and upper boundary
BC_node2_u=node_indx_u(mod(node_indx_u,(2*nx-1))==0|mod(node_indx_u,(2*nx-1))==1);  % left and right boundary
bound_u=sort(union(BC_node1_u, BC_node2_u));

Snunk=2*(2*nx-1)*(2*ny-1)+nx*ny;
for i_node=1:size(bound_u, 2);
    node=bound_u(1, i_node);
    ju=node;
    jv=node+(2*nx-1)*(2*ny-1);     
    
    x=node_xy(1,node);    y=node_xy(2,node);
    f(1:Snunk,1)=f(1:Snunk, 1)-a_S_linear(1:Snunk,ju).*True_solution(1,0,0,x,y,t);     
    a_S_linear(:,ju)=0.0E+00;      a_S_linear(ju,:)=0.0E+00;         
    
    a_S_linear(ju,ju)=1.0E+00;
    f(ju)=True_solution(1,0,0,x,y,t);  
    
    %p_S_linear(:,ju)=0.0E+00;      p_S_linear(ju,:)=0.0E+00;         
    %p_S_linear(ju,ju)=1.0E+00;
  
    
    f(1:Snunk,1)=f(1:Snunk, 1)-a_S_linear(1:Snunk,jv).*True_solution(2,0,0,x,y,t);     
    a_S_linear(:,jv)=0.0E+00;     a_S_linear(jv,:)=0.0E+00;        

    a_S_linear(jv,jv)=1.0E+00;
    f(jv)=True_solution(2,0,0,x,y,t);
    
    %p_S_linear(:,jv)=0.0E+00;     p_S_linear(jv,:)=0.0E+00;        
    %p_S_linear(jv,jv)=1.0E+00;
end

a_S=a_S_linear;
%p_S=p_S_linear;
f_S=f;