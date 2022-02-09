function b=Stokes_vector(ph,current_time,M_partition,T_partition,T_basis_u,number_of_elements,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vector_size_u,vector_size_p) 
%jv20191031

global lambda
number_of_local_basis_u=6;
number_of_local_basis_p=3;
T_basis_p=T_partition;
b1=assemble_vector_from_volume_integral_time_triangle('function_fu',current_time,M_partition,T_partition,T_basis_u,number_of_local_basis_u,number_of_elements,vector_size_u,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,0);
b2=assemble_vector_from_volume_integral_time_triangle('function_fv',current_time,M_partition,T_partition,T_basis_u,number_of_local_basis_u,number_of_elements,vector_size_u,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,0);
%b3=assemble_vector_from_volume_integral_time_triangle('function_fg',current_time,M_partition,T_partition,T_basis_p,number_of_local_basis_p,number_of_elements,vector_size_p,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0);
b4=assemble_vector_from_volume_integral_FE_triangle('FE_solution_triangle',ph,M_partition,T_partition,T_basis_p,T_basis_p,number_of_local_basis_p,number_of_elements,vector_size_p,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1,0,0);

b=[b1;b2;-1/lambda*b4];