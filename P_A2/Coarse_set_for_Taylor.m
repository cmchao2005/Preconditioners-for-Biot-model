function Coarse_set_for_Taylor(left,right,bottom,top,h_partition,dt)
%jv 20191119

global Coarse_prob lambda mu alpha c0 KK

N1_partition=(right-left)/h_partition(1);
N2_partition=(top-bottom)/h_partition(2);

nx=N1_partition+1;
ny=N2_partition+1;


%Mesh information for partition and finite element basis functions.
[M_partition,T_partition]=generate_M_T_triangle(left,right,bottom,top,h_partition,1);
Coarse_prob.node_xy_kp=M_partition;
Coarse_prob.element_node_kp=T_partition;
T_basis_kp=T_partition;
[M_basis_u,T_basis_u]=generate_M_T_triangle(left,right,bottom,top,h_partition,2);
Coarse_prob.node_xy_u=M_basis_u;
Coarse_prob.element_node_u=T_basis_u;

number_of_FE_nodes_u=(2*nx-1)*(2*ny-1);
number_of_FE_nodes_kp=nx*ny;

%set the boundary
[boundary_nodes_u,boundary_edges_u]=generate_boundary_nodes_edges(2*N1_partition,2*N2_partition,N1_partition,N2_partition);
bound_u=boundary_nodes_u(2,boundary_nodes_u(1,:)==-1);

[boundary_nodes_p,boundary_edges_p]=generate_boundary_nodes_edges(N1_partition,N2_partition,N1_partition,N2_partition);
bound_p=boundary_nodes_p(2,boundary_nodes_p(1,:)==-1);

bound_uvp=[bound_u,bound_u+number_of_FE_nodes_u,bound_p+2*number_of_FE_nodes_u+number_of_FE_nodes_kp];

Coarse_prob.bound_u=bound_u;
Coarse_prob.bound_p=bound_p;

%node_unknowns_u=setdiff(node_indx_u,bound_u);
%Coarse_prob.unkonwn_u=node_unknowns_u;
Coarse_prob.number_of_node_kp=number_of_FE_nodes_kp;

number_of_local_basis_u=6;
number_of_local_basis_kp=3;

number_of_elements=2*N1_partition*N2_partition;
matrix_size_uu=[number_of_FE_nodes_u number_of_FE_nodes_u];
matrix_size_ukp=[number_of_FE_nodes_u number_of_FE_nodes_kp];
matrix_size_kpkp=[number_of_FE_nodes_kp number_of_FE_nodes_kp];

[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(9);
Aup=sparse(matrix_size_ukp(1),matrix_size_ukp(2));
Auu=sparse(matrix_size_uu(1),matrix_size_uu(2));
App=sparse(matrix_size_kpkp(1),matrix_size_kpkp(2));

A1=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,1,0,2,1,0);
A2=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,1,2,0,1);
A3=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,1,0,2,0,1);
A4=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,1,2,1,0);
A5=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_u,number_of_local_basis_kp,number_of_local_basis_u,number_of_elements,matrix_size_ukp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,2,1,0);
A6=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_u,number_of_local_basis_kp,number_of_local_basis_u,number_of_elements,matrix_size_ukp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,2,0,1);
A7=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_kp,number_of_local_basis_kp,number_of_local_basis_kp,number_of_elements,matrix_size_kpkp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1,0,0);
A8=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_kp,number_of_local_basis_kp,number_of_local_basis_kp,number_of_elements,matrix_size_kpkp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,1,0,1,1,0);
A9=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_kp,number_of_local_basis_kp,number_of_local_basis_kp,number_of_elements,matrix_size_kpkp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,1,1,0,1);


Coarse_At=[2*mu*A1+mu*A2  mu*A3  -A5 Aup;
           mu*A4  2*mu*A2+mu*A1  -A6 Aup;
           -A5'  -A6'  -1/lambda*A7 alpha/lambda*A7;
            Aup' Aup' alpha/lambda*A7 -((c0+(alpha*alpha)/lambda)*A7+dt*KK*(A8+A9))];

Coarse_At(bound_uvp,:)=[];   Coarse_At(:,bound_uvp)=[];  
Coarse_prob.At=Coarse_At;

% M_Biot=[Auu  Auu  Aup Aup;
% Auu  Auu  Aup Aup;
% Aup' Aup' App App;
% Aup' Aup' alpha/lambda*A7 -(c0+(alpha*alpha)/lambda)*A7];
% 
% element_area = area_set(M_basis_u, number_of_elements, T_basis_u );
% 
% ph=zeros(number_of_FE_nodes_kp,1);
% ksi=zeros(number_of_FE_nodes_kp,1);
% for ii=1:number_of_FE_nodes_kp
% ph(ii)=True_solution(4,0,0,M_partition(1,ii),M_partition(2,ii),0);
% ksi(ii)=True_solution(3,0,0,M_partition(1,ii),M_partition(2,ii),0);
% end
% 
% bu=zeros(number_of_FE_nodes_u,1);
% vector_size_kp=number_of_FE_nodes_kp;
% vector_size_u=number_of_FE_nodes_u;
% for timestep=1:1
% 
%     current_time=dt*timestep;
% 
%     b1=assemble_vector_from_volume_integral_time_triangle('function_fu',current_time,M_partition,T_partition,T_basis_u,number_of_local_basis_u,number_of_elements,vector_size_u,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,0);
%     b2=assemble_vector_from_volume_integral_time_triangle('function_fv',current_time,M_partition,T_partition,T_basis_u,number_of_local_basis_u,number_of_elements,vector_size_u,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,0);
%     b3=zeros(number_of_FE_nodes_kp,1);
%     b4=assemble_vector_from_volume_integral_time_triangle('function_fphi',current_time,M_partition,T_partition,T_basis_kp,number_of_local_basis_kp,number_of_elements,vector_size_kp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0);
%     
%     b2=treat_Neumann_boundary_triangle('function_h2_N',b2,boundary_edges_u,M_partition,T_partition,T_basis_u,number_of_local_basis_u,current_time,2,0,0);
%     b4=treat_Neumann_boundary_triangle('function_g_N',b4,boundary_edges_p,M_partition,T_partition,T_basis_kp,number_of_local_basis_kp,current_time,1,0,0);
%     
%     b_vec=[b1;b2;b3;dt*b4];
%     b_up=[bu;bu;ksi;ph];
%     b_Biot=b_vec+M_Biot*b_up;
% 
% 
% 
%     [Fine_At,Fine_bt]=bc_treat_Biot(nx,ny,M_basis_u,Coarse_At,M_partition,b_Biot,current_time,bound_u,bound_p);
%     xst(bound_uvp)=Fine_bt(bound_uvp);
%     Fine_At(bound_uvp,:)=[];Fine_At(:,bound_uvp)=[];
%     Fine_bt(bound_uvp,:)=[];
% 
%     Fine_prob.ft=Fine_bt;
%     %         tol=1e-6; 
%     %         ini_sol=zeros(size(Fine_bt));
%     %         [st,flag,relres,iter,resvec]=gmres('AtimesX',Fine_bt, 400, tol, 500,'Pre_timesX',[], ini_sol);
%     %     relres_count(i_loop)=relres;
%     %     iter_count1(i_loop)=iter(1,2);
%     %     iter_count2(i_loop)=iter(1,1);
%         st=Fine_At\Fine_bt;
%         uvpnodes=[1:2*number_of_FE_nodes_u+2*number_of_FE_nodes_kp];
%         uvp_unkown=setdiff(uvpnodes,bound_uvp);
%         xst(uvp_unkown)=st;
%         ksi(:,1)=xst(1+2*number_of_FE_nodes_u:2*number_of_FE_nodes_u+number_of_FE_nodes_kp);
%         ph(:,1)=xst(1+2*number_of_FE_nodes_u+number_of_FE_nodes_kp:end);
%     fprintf(1,'   The timestep is %d \n',timestep);
% end
% 
%     uh1=xst(1:number_of_FE_nodes_u); 
%     [e_ul2,e_uh1]=errors_SU_new(element_area,T_basis_u,M_basis_u, uh1, number_of_elements, 6, 6,current_time); 
%     uh2=xst(1+number_of_FE_nodes_u:2*number_of_FE_nodes_u);
%     [e_vl2,e_vh1]=errors_SV_new(element_area,T_basis_u,M_basis_u, uh2, number_of_elements, 6, 6,current_time);
%     ksi=xst(1+2*number_of_FE_nodes_u:end);
%     [e_ksil2,e_ksih1]=errors_Sksi_new(element_area,T_partition,M_partition, ksi, number_of_elements, 3, 6,current_time);
%     [e_pl2,e_ph1]=errors_SP_new(element_area,T_partition,M_partition, ph, number_of_elements, 3, 6,current_time);
%     %fprintf(1, '# FineNX =  :');       fprintf(1, '  %d \n', FineNx);   % which method we are using?
%     fprintf(1, '# e_ul2 =  :');        fprintf(1, '  %d \n', e_ul2);
%     fprintf(1, '# e_uh1 =  :');        fprintf(1, '  %d \n', e_uh1);
%    fprintf(1, '# e_vl2 =  :');        fprintf(1, '  %d \n', e_vl2);
%     fprintf(1, '# e_vh1 =  :');        fprintf(1, '  %d \n', e_vh1);
%     fprintf(1, '# e_ksil2 =  :');      fprintf(1, '  %d \n', e_ksil2);
%     fprintf(1, '# e_ksih1 =  :');      fprintf(1, '  %d \n', e_ksih1);
%     fprintf(1, '# e_pl2 =  :');        fprintf(1, '  %d \n', e_pl2);
%     fprintf(1, '# e_ph1 =  :');        fprintf(1, '  %d \n', e_ph1);



