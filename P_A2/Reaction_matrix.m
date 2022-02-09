function [A, M]=Reaction_matrix(dt,M_partition,T_partition,number_of_elements,number_of_FE_nodes_p, Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle)
%j20191031
global KK c0 lambda
T_basis_p=T_partition;
matrix_size_pp=[number_of_FE_nodes_p,number_of_FE_nodes_p];
number_of_local_basis_p=3;
A1=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_p,T_basis_p,number_of_local_basis_p,number_of_local_basis_p,number_of_elements,matrix_size_pp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1,0,0);
A2=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_p,T_basis_p,number_of_local_basis_p,number_of_local_basis_p,number_of_elements,matrix_size_pp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,1,0,1,1,0);
A3=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_p,T_basis_p,number_of_local_basis_p,number_of_local_basis_p,number_of_elements,matrix_size_pp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,1,1,0,1);


A=(c0+1/lambda)*A1+dt*KK*(A2+A3);
M=(c0+1/lambda)*A1;