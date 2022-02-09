function A=Stokes_matrix(M_partition,T_partition,M_basis_u,T_basis_u,number_of_elements,number_of_FE_nodes_u,number_of_FE_nodes_ksi,...
       Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle)
%JV 20191031 

global mu lambda
   
matrix_size_uu=[number_of_FE_nodes_u number_of_FE_nodes_u];
matrix_size_up=[number_of_FE_nodes_u number_of_FE_nodes_ksi];
matrix_size_pp=[number_of_FE_nodes_ksi number_of_FE_nodes_ksi];
number_of_local_basis_u=6;
number_of_local_basis_p=3;
T_basis_p=T_partition;

A1=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,1,0,2,1,0);
A2=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,1,2,0,1);
A3=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,1,0,2,0,1);
A4=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,1,2,1,0);
A5=assemble_matrix_from_volume_integral_triangle('function_negativeone',M_partition,T_partition,T_basis_p,T_basis_u,number_of_local_basis_p,number_of_local_basis_u,number_of_elements,matrix_size_up,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,2,1,0);
A6=assemble_matrix_from_volume_integral_triangle('function_negativeone',M_partition,T_partition,T_basis_p,T_basis_u,number_of_local_basis_p,number_of_local_basis_u,number_of_elements,matrix_size_up,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,2,0,1);
A7=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_p,T_basis_p,number_of_local_basis_p,number_of_local_basis_p,number_of_elements,matrix_size_pp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1,0,0);

A=[2*mu*A1+mu*A2  mu*A3  A5;
   mu*A4  2*mu*A2+mu*A1  A6;
   A5'  A6'  -1/lambda*A7 ];


