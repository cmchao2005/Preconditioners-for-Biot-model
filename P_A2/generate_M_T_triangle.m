function [M,T]=generate_M_T_triangle(xl,xr,yb,yt,h_partition,basis_type)
%JvGuolaing, 07/10/2019.

nx=(xr-xl)/h_partition(1)+1;
ny=(yt-yb)/h_partition(2)+1;





if basis_type==2  
   M=zeros(2,(2*nx-1)*(2*ny-1));
   T=zeros(6,(nx-1)*(ny-1));
    for j = 1 : 2*ny-1
      for i = 1 : 2*nx - 1
           M(1,i+(j-1)*(2*nx-1)) =    ...
                ( ( 2 * nx - i - 1 ) * xl   ...
                + (          i - 1 ) * xr ) ...
                / ( 2 * nx     - 2 );

            M(2,i+(j-1)*(2*nx-1)) =    ...
                ( ( 2 * ny - j - 1 ) * yb   ...
                + (          j - 1 ) * yt ) ...
                / ( 2 * ny     - 2 );

      end
    end
    
    element = 0;
    for j = 1 : ny - 1
      for i = 1 : nx - 1

          sw = ( j - 1 ) * 2 * ( 2 * nx - 1 ) + 2 * i - 1;
          s  = sw + 1;
          se = sw + 2;

          w  = sw + 2 * nx - 1;
          c  = w  + 1;
          e  = w  + 2;

          nw = w  + 2 * nx - 1;
          n  = nw + 1;
          ne = nw + 2;

          element = element + 1;
          T(1,element) = sw;
          T(2,element) = se;
          T(3,element) = nw;
          T(4,element) = s;
          T(5,element) = c;
          T(6,element) = w;

          element = element + 1;
          T(1,element) = ne;
          T(2,element) = nw;
          T(3,element) = se;
          T(4,element) = n;
          T(5,element) = c;
          T(6,element) = e;

      end
    end
    
else 
    M=zeros(2,nx*ny);
    T=zeros(3,(nx-1)*(ny-1));
    for j = 1 : ny
      for i = 1 : nx
           M(1,i+(j-1)*nx) =    ...
                ( ( nx - i ) * xl   ...
                + (  i - 1 ) * xr ) ...
                / (  nx  - 1 );

            M(2,i+(j-1)*nx) =    ...
                ( ( ny - j ) * yb   ...
                + ( j - 1 ) * yt ) ...
                / ( ny  - 1 );

      end
     end
   
    element = 0;
    for j = 1 : ny - 1
      for i = 1 : nx - 1

          sw = ( j - 1 ) *  nx  +  i ;
          se = sw + 1;
          nw = sw  + nx;
          ne = nw + 1;

          element = element + 1;
          T(1,element) = sw;
          T(2,element) = se;
          T(3,element) = nw;




          element = element + 1;
          T(1,element) = ne;
          T(2,element) = nw;
          T(3,element) = se;



       end
    end
end