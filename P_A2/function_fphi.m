function r=function_fphi(t,x,y)
%jv 20190930
global c0 KK alpha
r=-((-c0+2*KK)*sin(x+y)*exp(-t)-alpha*(cos(x)+cos(y))*exp(-t));