function element_area = area_set (node_xy, element_num, element_node )

%% AREA_SET sets the area of each element.(included in FEM_STOKES_DARCY)
%    The areas of the elements are needed in order to adjust
%    the integral estimates produced by the quadrature formulas.
%  Modified:
%    19 Jun 2014
%
%  Parameters:
%    Input:
%        integer NODE_NUM, the number of nodes.
%        real NODE_XY(2,NODE_NUM), the nodes.
%        integer NNODES, the number of local nodes per element.
%        integer ELEMENT_NUM, the number of elements.
%        integer ELEMENT_NODE(NNODES,ELEMENT_NUM);
%        ELEMENT_NODE(I,J) is the global index of local node I in element J.
%
%    Output:
%        real ELEMENT_AREA(ELEMENT_NUM), the area of elements.
%
  for element = 1 : element_num

    i1 = element_node(1,element);
    x1 = node_xy(1,i1);
    y1 = node_xy(2,i1);

    i2 = element_node(2,element);
    x2 = node_xy(1,i2);
    y2 = node_xy(2,i2);

    i3 = element_node(3,element);
    x3 = node_xy(1,i3);
    y3 = node_xy(2,i3);

    element_area(element) = 0.5E+00 * abs ...
      ( y1 * ( x2 - x3 ) ...
      + y2 * ( x3 - x1 ) ...
      + y3 * ( x1 - x2 ) );

  end
