function main_Hh
%jv 20191119

clc;
clear;
global lambda mu Coarse_prob Fine_prob enforce_local_0Diri Nsubdomains level c0 alpha KK detau
init_time=0;
end_time=1e-5;
dt=1e-5;

left=0;right=1;bottom=0;top=1;

c0=1; KK=1;  alpha=1; level=2; enforce_local_0Diri=1;
E=1;nu=0.3;
%lambda =1000 ; mu=1;
lambda = E*nu/((1+nu)*(1-2*nu)); mu= E/(2*(1+nu));
%hpart=8;
log_filename=strcat('N',mat2str(level),'_Hh_',mat2str(nu),'_run_log', '.txt');    
fid = fopen(log_filename, 'wt'); 

fprintf(1,'KK=%d, level=%d, lambda=%d, \n',[KK,level,lambda]);
%for ii=1:9
for detau=2:6
for npart=2:2
     %detau=64/(2^ii)
     hpart=8*detau;
    Coarse_prob={};
    Fine_prob={};
    
    % basic set
    CoarseNx=npart; CoarseNy=CoarseNx;
    Nsubdomains=CoarseNx*CoarseNy;
    fprintf(1,'Nsubdomains=%d,  overlap=%d, \n',[Nsubdomains,detau]);
    fprintf(fid,'Nsubdomains=%d, overlap=%d, E=%f, nu=%f, level=%d, \n',[Nsubdomains,detau,E,nu,level]);
    FineNx=(hpart*npart); FineNy=FineNx;
    
    Coarse_prob.hx = 1/CoarseNx; Coarse_prob.hy = 1/CoarseNy;
    Fine_prob.hx=1/FineNx; Fine_prob.hy=1/FineNy;

    Coarse_set_for_Taylor(left,right,bottom,top,[Coarse_prob.hx,Coarse_prob.hy],dt)
    
    Biot_solver(fid,left,right,bottom,top,FineNx,init_time,end_time,dt);

    overlap_size=detau*Fine_prob.hx;
    Generate_Restriction(overlap_size, CoarseNx, CoarseNy, FineNx, FineNy,left, bottom, right, top);
    
    tol=1e-6; 
    ini_sol=zeros(size(Fine_prob.ft));
    
    [xst,flag,relres,iter,resvec]=gmres('AtimesX',Fine_prob.ft, 400, tol, 500,'PreDtimesX',[], ini_sol);
    fprintf(1,'   Restart=%d,  the iter=%d, and the relres =%e,\n\n',[iter(1),iter(2),relres]);
    fprintf(fid,'   Restart=%d, the iter=%d, and the relres =%e,\n\n',[iter(1),iter(2),relres]);
%     filea=['resvec1_c.txt'];
%     dlmwrite(filea,resvec,'delimiter','\t');
   
end
end
fprintf(1,'PA, \n');

% lambda=lambda*10
% end


