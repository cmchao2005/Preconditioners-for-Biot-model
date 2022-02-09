function r=function_h1_N(t,x,y)
%jv 20190930
global mu lambda alpha
r=2*mu*cos(x)*exp(-t)+lambda*(cos(x)+cos(y))*exp(-t)-alpha*sin(x+y)*exp(-t);