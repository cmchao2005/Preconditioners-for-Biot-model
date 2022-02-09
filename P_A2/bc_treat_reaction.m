function [a_S, f_S] = bc_treat_reaction(nx,ny,node_xy,a_S_linear,f,t)


node_indx_p=[1:nx*ny];
BC_node1_p=union([1:nx], [nx*ny-nx+1:nx*ny]);   % lower and upper boundary
BC_node2_p=node_indx_p(mod(node_indx_p,nx)==0|mod(node_indx_p,nx)==1);  % left and right boundary
bound_p=sort(union(BC_node1_p, BC_node2_p));

Snunk=nx*ny;
for i_node=1:size(bound_p, 2);
    node=bound_p(1, i_node);
    ju=node;  

    x=node_xy(1,node);    y=node_xy(2,node);
    f(1:Snunk,1)=f(1:Snunk, 1)-a_S_linear(1:Snunk,ju).*True_solution(4,0,0,x,y,t);     
    a_S_linear(:,ju)=0.0E+00;      a_S_linear(ju,:)=0.0E+00;         
    
    a_S_linear(ju,ju)=1.0E+00;
    f(ju)=True_solution(4,0,0,x,y,t);  
    
    
end

a_S=a_S_linear;
f_S=f;