function r=function_fv(t,x,y)
%jv 20190930
global mu lambda alpha
r=(lambda+2*mu)*sin(y)*exp(-t)+alpha*cos(x+y)*exp(-t);