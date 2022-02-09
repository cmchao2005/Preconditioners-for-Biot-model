function r=Gauss_quadrature_for_volume_integral_test_FE_triangle(coefficient_function_name,ph_local_1,Gauss_coefficient_local,Gauss_point_local,vertices,test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y,FE_basis_type,FE_derivative_degree_x,FE_derivative_degree_y)
%Xiaoming He, 02/07/2010.
%Use Gauss quadrature to compute a volume integral on a local triangular element T for a vector.
%We will use "FE" to replace "finite element" in the comments.
%The integrand of the volume integral must be in the following format:
%a coefficient function which consists of two FE solutions * a test FE basis function (or its derivatives).
%coefficient_function_name: the coefficient function of the integrand.
%Gauss_coefficient_local,Gauss_point_local:the Gauss coefficients and Gauss points on the triangular element T.
%vertices: the coordinates of all vertices of the triangular element T.
%test_basis_type:the type of the test FE basis function.
%test_basis_type=1:2D linear FE.  
%test_basis_type=2:2D Lagrange quadratic FE.
%test_basis_index: the index of test FE basis function to specify which test FE basis function we want to use.
%test_derivative_degree_x:the derivative degree of the test FE basis function with respect to x.
%test_derivative_degree_y:the derivative degree of the test FE basis function with respect to y.
%FE_basis_type:the type of the FE basis function used for the coefficient function.
%FE_basis_type=1:2D linear FE.  
%FE_basis_type=2:2D Lagrange quadratic FE.
%FE_basis_index: the index of the FE basis function used for the coefficient function to specify which test FE basis function we want to use.
%FE_derivative_degree_x:the derivative degree of the FE basis function used for the coefficient function with respect to x.
%FE_derivative_degree_y:the derivative degree of the FE basis function used for the coefficient function with respect to y.
%uh_local: the values of the FE solution at the nodes of FE in a triangular element.
%uh_local_1: the uh_local for the first FE solution in the coefficient function.
%uh_local_2: the uh_local for the second FE solution in the coefficient function.

%Gpn: the number of the Gauss points of the Gauss quadrature we are using.

Gpn=length(Gauss_coefficient_local);
r=0;
for i=1:Gpn
      r=r+Gauss_coefficient_local(i)...
          *feval(coefficient_function_name,Gauss_point_local(i,1),Gauss_point_local(i,2),ph_local_1,vertices,FE_basis_type,FE_derivative_degree_x,FE_derivative_degree_y)...
          *triangular_local_basis(Gauss_point_local(i,1),Gauss_point_local(i,2),vertices,test_basis_type,test_basis_index,test_derivative_degree_x,test_derivative_degree_y);
end