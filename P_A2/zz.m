clc;
clear;

global lambda mu  Fine_prob c0 alpha KK drop_tol

c0=1; KK=1;  alpha=1; 
E=1;nu=0.499;
%lambda =1000 ; mu=1;
lambda = E*nu/((1+nu)*(1-2*nu)); mu= E/(2*(1+nu));

fid=0;
init_time=0;
end_time=1e-5;   % we only computer 1 step
dt=1e-5;
FineNx=32;
drop_tol=1e-3;


left=0;right=1;bottom=0;top=1;

Biot_solver(fid,left,right,bottom,top,FineNx,init_time,end_time,dt)