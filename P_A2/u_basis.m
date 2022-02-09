function b=u_basis(j,r,s)
    if ( j == 1 )
     b    =   2.0E+00 *     ( 1.0E+00 - r - s ) * ( 0.5E+00 - r - s );
    elseif ( j == 2 )
     b    =   2.0E+00 * r * ( r - 0.5E+00 );   
    elseif ( j == 3 )
     b    =   2.0E+00 * s * ( s - 0.5E+00 );
    elseif ( j == 4 )
     b    =   4.0E+00 * r * ( 1.0E+00 - r - s );
    elseif ( j == 5 )
     b    =   4.0E+00 * r * s;
    elseif ( j == 6 )
     b    =   4.0E+00 * s * ( 1.0E+00 - r - s );
    end