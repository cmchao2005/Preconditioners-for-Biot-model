function [ el2, eh1 ] =  errors_SU_new ( element_area, element_node, node_xy, numerical_sol, element_num, nnodes, nq,t); 

%% ERRORS calculates the L2 and H1-seminorm errors.
%
%  Modified:
%
%    17 May 2005
%  Parameters:
%
%    Input, real ELEMENT_AREA(ELEMENT_NUM), the area of elements.
%
%    Input, integer ELEMENT_NODE(NNODES,ELEMENT_NUM);
%    ELEMENT_NODE(I,J) is the global index of local node I in element J.
%
%    Input, integer INDX(NODE_NUM), gives the index of the unknown quantity
%    associated with the given node.
%
%    Input, real NODE_XY(2,NODE_NUM), the nodes.
%
%    Input, real F(NUNK), the coefficients of the solution.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer NNODES, the number of nodes used to form one element.
%
%    Input, integer nq, the number of points in the quadrature rule.
%    This is actually fixed at 13.
%
%    Input, integer NUNK, the number of unknowns.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Output, real EL2, the L2 error.
%
%    Output, real EH1, the H1 seminorm error.
%
%  Local Parameters:
%
%    Local, double precision AR, the weight for a given quadrature point
%    in a given element.
%
%    Local, double precision BB, BX, BY, a basis function and its first
%    derivatives evaluated at a particular quadrature point.
%
%    Local, double precision EH1, the H1 seminorm error.
%
%    Local, double precision EL2, the L2 error.
%
%    Local, double precision UEX, UEXX, UEXY, the exact solution and its first
%    derivatives evaluated at a particular quadrature point.
%
%    Local, double precision UH, UHX, UHY, the computed solution and its first
%    derivatives evaluated at a particular quadrature point.
%
%    Local, double precision WQE(nq), stores the quadrature weights.
%
%    Local, double precision X, Y, the coordinates of a particular
%    quadrature point.
%
%    Local, double precision XQE(nq), YQE(nq), stores the location
%    of quadrature points in a given element.
%

  el2 = 0.0E+00;
  eh1 = 0.0E+00;
  
%  ul2 = 0.0E+00;
%  uh1 = 0.0E+00;
  
  [ wq, quad_xy ] = quad_rule (nq);     % For reference element
  for element = 1:element_num

    t3(1:2,1:3) = node_xy(1:2,element_node(1:3,element));     % 4 nodes in element
%
%  Map the quadrature points QUAD_XY to points XY in the physical triangle.
%
    xy = reference_to_physical_t3( t3, nq, quad_xy );    

    for quad = 1 : nq

      x = xy(1,quad);
      y = xy(2,quad);
      w = element_area(element) * wq(quad);
      
      uh = 0.0E+00;
      dudxh = 0.0E+00;
      dudyh = 0.0E+00;

      for node = 1 : nnodes
          
          iu = element_node(node,element);
%%          ip=iu-(2*nx-1)^2;
          
          
          [ bi, dbidx, dbidy ] = qbf ( x, y, element, node, node_xy, element_node);
%        i = indx(ip);    % node_u_variable, node_v_variable,
%        node_p_variable.
          uh    = uh    + bi    * numerical_sol(iu);
          dudxh = dudxh + dbidx * numerical_sol(iu);
          dudyh = dudyh + dbidy * numerical_sol(iu);
          
      end
      u = True_solution(1,0,0,x,y,t);    % Provide exact solution at (xq,yq)
      dudx=True_solution(1,1,0,x,y,t);
      dudy =True_solution(1,0,1,x,y,t);
      el2 = el2 + ( uh - u )^2 * w;                                  % Second order derivatives
      eh1 = eh1 + ( ( dudxh - dudx )^2 + ( dudyh - dudy )^2 ) * w;
      
%      ul2 = ul2 + u^2*w;
%      uh1 = uh1 + (dudx^2+dudy^2)*w;
      
    end
  end
  
  el2  = sqrt ( el2 );
  eh1  = sqrt ( eh1 );
   
  %rel2 = sqrt ( el2 )/sqrt(ul2);    %%%    Relative error
  %reh1 = sqrt ( eh1 )/sqrt(uh1);    %%%    Relative error

  fprintf ( 1, '\n' );
  fprintf ( 1, '*********************************************\n' );
  fprintf ( 1, '*                                           *\n' );
  fprintf ( 1, '*  ERRORS_u:                              *\n' );
  fprintf ( 1, '*    L2 error =          %e       *\n', el2 );
  fprintf ( 1, '*    H1-seminorm error = %e       *\n', eh1 );
%  fprintf ( 1, '*  Relative L2 error =   %16.15f      *\n', rel2 );
%  fprintf ( 1, '*  Relative H1-seminorm error = %16.15f   *\n', reh1 );
  fprintf ( 1, '*                                           *\n' );
  fprintf ( 1, '*********************************************\n' );
  