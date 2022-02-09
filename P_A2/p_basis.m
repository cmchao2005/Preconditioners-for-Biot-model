function bp=p_basis(j,r,s)
    if ( j == 1 )
        bp  = (1.0E+00 - r - s );
    elseif ( j == 2 )
        bp  = r ;
    elseif ( j == 3 )
        bp  = s; 
    end