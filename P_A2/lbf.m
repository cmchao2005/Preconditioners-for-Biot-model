function [ b, dbdx, dbdy ] = lbf ( x, y, element, inode, node_xy, element_node)

%% LBF evaluates the linear basis functions.
%
%  Discussion:
%
%    This routine assumes that the "midpoint" nodes are, in fact,
%    exactly the average of the two extreme nodes.  
%
%    Assuming this property of the midpoint nodes makes it easy to
%    determine the values of (R,S) in the reference element that
%    correspond to (X,Y) in the physical element.
%
%    Once we know the (R,S) coordinates, it's easy to evaluate the
%    basis functions and derivatives.
%
%  The physical element T3:
%
%    In this picture, we don't mean to suggest that the bottom of
%    the physical triangle is horizontal.  However, we do assume that
%    each of the sides is a straight line, and that the intermediate
%    points are exactly halfway on each side.
%
%    |
%    |
%    |        3
%    |       / \ 
%    |      /   \
%    |     /     \
%    |    1-------2
%    |
%    +--------X-------->
%
%  Reference element T6:
%
%    In this picture of the reference element, we really do assume
%    that one side is vertical, one horizontal, of length 1.
%
%    |
%    |
%    1  3
%    |  |\
%    |  | \
%    |  |  \
%    |  |   \
%    0  1----2
%    |
%    +--0--R--1-------->
%
%  Modified:
%
%    09 Sep 2005
%
%  Parameters:
%
%    Input, real X, Y, the (global) coordinates of the point
%    at which the basis function is to be evaluated.
%
%    Input, integer ELEMENT, the index of the element which contains the point.
%
%    Input, integer INODE, the local index (between 1 and 6) that
%    specifies which basis function is to be evaluated.
%
%    Input, real NODE_XY(2,NODE_NUM), the nodes.
%
%    Input, integer ELEMENT_NODE(NNODES,ELEMENT_NUM);
%    ELEMENT_NODE(I,J) is the global index of local node I in element J.
%
%    Input, integer ELEMENT_NUM, the number of elements.
%
%    Input, integer NNODES, the number of nodes used to form one element.
%
%    Input, integer NODE_NUM, the number of nodes.
%
%    Output, real B, DBDX, DBDY, the value of the basis function
%    and its X and Y derivatives at (X,Y).
% 
%  The coordinates of 6 nodes in element
  xn(1:3) = node_xy(1,element_node(1:3,element));
  yn(1:3) = node_xy(2,element_node(1:3,element));
%
%  Determine the (R,S) coordinates corresponding to (X,Y).
%
%  What is happening here is that we are solving the linear system:
%
%    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
%    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
%
%  by computing the inverse of the coefficient matrix and multiplying
%  it by the right hand side to get R and S.
%
%  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
%  for R and S.
%
  det =   ( xn(2) - xn(1) ) * ( yn(3) - yn(1) ) ...
        - ( xn(3) - xn(1) ) * ( yn(2) - yn(1) );

  r = ( ( yn(3) - yn(1) ) * ( x     - xn(1) ) ...
      + ( xn(1) - xn(3) ) * ( y     - yn(1) ) ) / det;

  drdx = ( yn(3) - yn(1) ) / det;
  drdy = ( xn(1) - xn(3) ) / det;

  s = ( ( yn(1) - yn(2) ) * ( x     - xn(1) ) ...
      + ( xn(2) - xn(1) ) * ( y     - yn(1) ) ) / det;

  dsdx = ( yn(1) - yn(2) ) / det;
  dsdy = ( xn(2) - xn(1) ) / det;
%
%  The basis functions can now be evaluated in terms of the
%  reference coordinates R and S.  It's also easy to determine
%  the values of the derivatives with respect to R and S.
%
  if ( inode == 1 )

    b    =   (1.0E+00 - r - s );
    dbdr = - 1.0E+00;
    dbds = - 1.0E+00;

  elseif ( inode == 2 )

    b    =   r;
    dbdr =   1.0E+00;
    dbds =   0.0E+00;

  elseif ( inode == 3 )

    b    =   s ;
    dbdr =   0.0E+00;
    dbds =   1.0E+00;
  else

    fprintf ( 1, '\n' );
    fprintf ( 1, 'LBF - Fatal error!\n' );
    fprintf ( 1, '  Request for local basis function INODE = %d\n', inode );
    error ( 'LBF - Fatal error!' );

  end
%
%  We need to convert the derivative information from (R(X,Y),S(X,Y))
%  to (X,Y) using the chain rule.
%
  dbdx = dbdr * drdx + dbds * dsdx;
  dbdy = dbdr * drdy + dbds * dsdy;
