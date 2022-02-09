function [a_S, f_S] =bc_treat_Biot(nx,ny,node_xy_u,a_S_linear,node_xy_p,f,t,bound_u,bound_p)
%jv 20191119



Snunk=2*(2*nx-1)*(2*ny-1)+2*nx*ny;
for i_node=1:size(bound_u, 2);
    node=bound_u(1, i_node);
    ju=node;
    jv=node+(2*nx-1)*(2*ny-1);     
    
    x=node_xy_u(1,node);    y=node_xy_u(2,node);
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

for i_node_p=1:size(bound_p, 2);
    node_p=bound_p(1, i_node_p);
    jp=node_p+2*(2*nx-1)*(2*ny-1)+nx*ny;    
    
    x=node_xy_p(1,node_p);    y=node_xy_p(2,node_p);
    f(1:Snunk,1)=f(1:Snunk, 1)-a_S_linear(1:Snunk,jp).*True_solution(4,0,0,x,y,t);     
    a_S_linear(:,jp)=0.0E+00;      a_S_linear(jp,:)=0.0E+00;         
    
    a_S_linear(jp,jp)=1.0E+00;
    f(jp)=True_solution(4,0,0,x,y,t);  
end

a_S=a_S_linear;
f_S=f;