function b=Reaction_vector(current_time,dt,ph,dksi,M_Reaction,M_partition,T_partition,number_of_elements,vector_size_p, Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle)

global lambda

T_basis_p=T_partition;
number_of_local_basis_p=3;
%b1=assemble_vector_from_volume_integral_FE_triangle('FE_solution_triangle',ph,M_partition,T_partition,T_basis_p,T_basis_p,number_of_local_basis_p,number_of_elements,vector_size_p,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1,0,0);
b2=assemble_vector_from_volume_integral_FE_triangle('FE_solution_triangle',dksi,M_partition,T_partition,T_basis_p,T_basis_p,number_of_local_basis_p,number_of_elements,vector_size_p,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1,0,0);
b3=assemble_vector_from_volume_integral_time_triangle('function_fphi',current_time,M_partition,T_partition,T_basis_p,number_of_local_basis_p,number_of_elements,vector_size_p,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0);

b=dt*b3+1/lambda*b2+M_Reaction*ph;