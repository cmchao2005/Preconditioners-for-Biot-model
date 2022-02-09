function r=function_fu(t,x,y)
%jv 20190930
global mu lambda alpha
r=(lambda+2*mu)*sin(x)*exp(-t)+alpha*cos(x+y)*exp(-t);