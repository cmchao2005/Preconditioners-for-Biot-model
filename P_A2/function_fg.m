function r=function_fg(t,x,y)
%jv 20190930
global mu lambda
r=-1/lambda*sin(x+y)*exp(-t);