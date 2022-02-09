function Biot_solver(fid,left,right,bottom,top,FineNx,init_time,end_time,dt)
% jv 20191031

global mu lambda c0 KK Fine_prob alpha chol_Au chol_SSksi chol_SSp drop_tol

Nt=(end_time-init_time)/dt;
FineNy=FineNx;

h_partition=[1/FineNx,1/FineNy];

[M_partition,T_partition]=generate_M_T_triangle(left,right,bottom,top,h_partition,1);
[M_basis_u,T_basis_u]=generate_M_T_triangle(left,right,bottom,top,h_partition,2);
Fine_prob.node_xy_kp=M_partition;     %kp means ksi and P which have same structure
Fine_prob.element_node_kp=T_partition;
Fine_prob.node_xy_u=M_basis_u;
Fine_prob.element_node_u=T_basis_u;

number_of_elements=2*FineNx*FineNy;
number_of_FE_nodes_u=(2*FineNx+1)*(2*FineNy+1);
number_of_FE_nodes_kp=(FineNx+1)*(FineNy+1);


[boundary_nodes_u,boundary_edges_u]=generate_boundary_nodes_edges(2*FineNx,2*FineNy,FineNx,FineNy);
bound_u=boundary_nodes_u(2,boundary_nodes_u(1,:)==-1);
Numan_bound_u=boundary_nodes_u(2,boundary_nodes_u(1,:)==-2);
Fine_prob.bound_u=bound_u;
Fine_prob.Numan_bound_u=Numan_bound_u;


[boundary_nodes_p,boundary_edges_p]=generate_boundary_nodes_edges(FineNx,FineNy,FineNx,FineNy);
bound_p=boundary_nodes_p(2,boundary_nodes_p(1,:)==-1);
Numan_bound_p=boundary_nodes_p(2,boundary_nodes_p(1,:)==-2);
Fine_prob.bound_p=bound_p;
Fine_prob.Numan_bound_p=Numan_bound_p;

% u_value(1:number_of_FE_nodes_u)=0;
% p_value(1:number_of_FE_nodes_kp)=0;
% u_value(bound_u)=1;u_value(Numan_bound_u)=0.5;
% p_value(bound_p)=1;p_value(Numan_bound_p)=0.5;
% BBBu=reshape(u_value',[2*FineNx+1,2*FineNx+1]);
% BBBp=reshape(p_value',[FineNx+1,FineNx+1]);
% BBBu'
% BBBp'

bound_uvp=[bound_u,bound_u+number_of_FE_nodes_u,bound_p+2*number_of_FE_nodes_u+number_of_FE_nodes_kp];

Fine_prob.number_of_unknown_u=number_of_FE_nodes_u-size(bound_u,2);
Fine_prob.number_of_unknown_p=number_of_FE_nodes_kp-size(bound_p,2);
element_area = area_set(M_basis_u, number_of_elements, T_basis_u );
[Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle]=generate_Gauss_reference_triangle(9);
matrix_size_uu=[number_of_FE_nodes_u number_of_FE_nodes_u];
matrix_size_ukp=[number_of_FE_nodes_u number_of_FE_nodes_kp];
matrix_size_kpkp=[number_of_FE_nodes_kp number_of_FE_nodes_kp];
number_of_local_basis_u=6;
number_of_local_basis_kp=3;
T_basis_kp=T_partition;
vector_size_kp=number_of_FE_nodes_kp;
vector_size_u=number_of_FE_nodes_u;

Fine_prob.number_of_node_u=number_of_FE_nodes_u;
Fine_prob.number_of_node_kp=number_of_FE_nodes_kp;

Aup=sparse(matrix_size_ukp(1),matrix_size_ukp(2));
Auu=sparse(matrix_size_uu(1),matrix_size_uu(2));
App=sparse(matrix_size_kpkp(1),matrix_size_kpkp(2));
bu=zeros(number_of_FE_nodes_u,1);

A1=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,1,0,2,1,0);
A2=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,1,2,0,1);
A3=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,1,0,2,0,1);
A4=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_u,T_basis_u,number_of_local_basis_u,number_of_local_basis_u,number_of_elements,matrix_size_uu,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,1,2,1,0);
A5=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_u,number_of_local_basis_kp,number_of_local_basis_u,number_of_elements,matrix_size_ukp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,2,1,0);
A6=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_u,number_of_local_basis_kp,number_of_local_basis_u,number_of_elements,matrix_size_ukp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,2,0,1);
A7=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_kp,number_of_local_basis_kp,number_of_local_basis_kp,number_of_elements,matrix_size_kpkp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0,1,0,0);
A8=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_kp,number_of_local_basis_kp,number_of_local_basis_kp,number_of_elements,matrix_size_kpkp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,1,0,1,1,0);
A9=assemble_matrix_from_volume_integral_triangle('function_one',M_partition,T_partition,T_basis_kp,T_basis_kp,number_of_local_basis_kp,number_of_local_basis_kp,number_of_elements,matrix_size_kpkp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,1,1,0,1);

A_Biot=[2*mu*A1+mu*A2  mu*A3  -A5 Aup;
        mu*A4  2*mu*A2+mu*A1  -A6 Aup;
        -A5'   -A6'  -1/lambda*A7 alpha/lambda*A7;
         Aup' Aup' alpha/lambda*A7 -((c0+(alpha*alpha)/lambda)*A7+dt*KK*(A8+A9))];
     
     

M_Biot=[Auu  Auu  Aup Aup;
Auu  Auu  Aup Aup;
Aup' Aup' App App;
Aup' Aup' alpha/lambda*A7 -(c0+(alpha*alpha)/lambda)*A7];

bound_uv= [bound_u,bound_u+number_of_FE_nodes_u];   
     
A_u=[2*mu*A1+mu*A2  mu*A3;
     mu*A4  2*mu*A2+mu*A1]; 
A_u(bound_uv,:)=[];A_u(:,bound_uv)=[];

B_uksi=[-A5', -A6'];
B_uksi(:,bound_uv)=[];

A_ksi=1/lambda*A7;
B_pksi=alpha/lambda*A7;
B_pksi(bound_p,:)=[];
A_p=-((c0+(alpha*alpha)/lambda)*A7+dt*KK*(A8+A9));
A_p(:,bound_p)=[];A_p(bound_p,:)=[];

% S_ksi=A_ksi+B_uksi*inv(A_u)*B_uksi';  % Exact Schur complement
% S_p=A_p+B_pksi*inv(S_ksi)*B_pksi';    % Exact Schur complement

SS_ksi = (2*mu+lambda)/(2*mu*lambda)*A7;    %MCai, SS_ksi is SPD, analytic app
%MCai, SS_p is negative PD, analytic app
SS_p = -((c0+(alpha*alpha)/lambda)*A7+dt*KK*(A8+A9))+(2*mu*alpha*alpha)/(lambda*(lambda+2*mu))*A7; 
SS_p(:,bound_p)=[];SS_p(bound_p,:)=[];

B_up=[Aup;Aup];
B_up(bound_uv,:)=[];
Bpp=App;
Bpp(:,bound_p)=[];

B_up3=B_up;
B_up3(:,bound_p)=[];

%三步收敛到真解
% P_Biot=[A_u,     B_up,  B_up3;
%         B_uksi, -S_ksi, Bpp;
%         B_up3',  B_pksi, S_p];
    
% Fine_prob.Pt= P_Biot;


    
Fine_prob.Au=A_u;
Fine_prob.SS_ksi=SS_ksi;  %正定
Fine_prob.SS_p = SS_p;    %负定
Fine_prob.B_uksi=B_uksi;
Fine_prob.B_pksi=B_pksi;

%理想状况
% P_Biot=[A_u,      B_uksi,    B_up3;
%         B_uksi',  -SS_ksi,   Bpp;
%         B_up3',  Bpp',      -SS_p];

ph=zeros(number_of_FE_nodes_kp,1);
ksi=zeros(number_of_FE_nodes_kp,1);
for ii=1:number_of_FE_nodes_kp
ph(ii)=True_solution(4,0,0,M_partition(1,ii),M_partition(2,ii),init_time);
ksi(ii)=True_solution(3,0,0,M_partition(1,ii),M_partition(2,ii),init_time);
end
 
for timestep=1:Nt

    current_time=init_time+dt*timestep;

    b1=assemble_vector_from_volume_integral_time_triangle('function_fu',current_time,M_partition,T_partition,T_basis_u,number_of_local_basis_u,number_of_elements,vector_size_u,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,0);
    b2=assemble_vector_from_volume_integral_time_triangle('function_fv',current_time,M_partition,T_partition,T_basis_u,number_of_local_basis_u,number_of_elements,vector_size_u,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,2,0,0);
    b3=zeros(number_of_FE_nodes_kp,1);
    b4=assemble_vector_from_volume_integral_time_triangle('function_fphi',current_time,M_partition,T_partition,T_basis_kp,number_of_local_basis_kp,number_of_elements,vector_size_kp,Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,1,0,0);
    b2=treat_Neumann_boundary_triangle('function_h2_N',b2,boundary_edges_u,M_partition,T_partition,T_basis_u,number_of_local_basis_u,current_time,2,0,0);
    b4=treat_Neumann_boundary_triangle('function_g_N',b4,boundary_edges_p,M_partition,T_partition,T_basis_kp,number_of_local_basis_kp,current_time,1,0,0);
    b_vec=[b1;b2;b3;dt*b4];
    b_up=[bu;bu;ksi;ph];
    b_Biot=b_vec+M_Biot*b_up;
    [Fine_At,Fine_bt]=bc_treat_Biot(FineNx+1,FineNy+1,M_basis_u,A_Biot,M_partition,b_Biot,current_time,bound_u,bound_p);
    xst(bound_uvp)=Fine_bt(bound_uvp);
    Fine_At(bound_uvp,:)=[];Fine_At(:,bound_uvp)=[];
    Fine_bt(bound_uvp,:)=[];
    Fine_prob.At= Fine_At;
 
	% added by MCai, for comparion, put a positive sign
	chol_Au=ichol(Fine_prob.Au, struct('type','ict','droptol',drop_tol));  
	chol_SSksi=ichol(Fine_prob.SS_ksi, struct('type','ict','droptol',drop_tol)); 
	chol_SSp= ichol(-Fine_prob.SS_p, struct('type','ict','droptol',drop_tol));  
    Fine_prob.ft=Fine_bt;
    tol=1e-6; 
    ini_sol=zeros(size(Fine_bt));
    
    fileID = fopen('results.txt','w');

%     fprintf(fileID, 'Block Diag with + sign in all blocks \n');    % P_D1
%     tStart = cputime;
%    [st,flag,relres,iter,resvec]=gmres('AtimesX',Fine_bt, 960, tol, 1200,'Pre_D1_ichol',[], ini_sol);
%      
%  	tEndD1 = cputime - tStart
%     flag 
%     iter
%     fprintf(fileID,'%5d %5d %6d\n', flag, iter);
%  
% %      Pre_D1=blkdiag(Fine_prob.Au, Fine_prob.SS_ksi, -Fine_prob.SS_p);    %correct, (2.4) 02/04
% %      [eigvs, D, flg]=eigs(inv(Pre_D1)*Fine_prob.At, 654);
% %      eig_vals=diag(D)
% %      flg
% %      fprintf(fileID,'%5d \n', eig_vals);
% %      figure
% %      plot(eig_vals, '^');    
% %      return
%     
%     fprintf(fileID,'Block Diag with - sign in S_p \n');                  % P_D2
%     tStart = cputime;
%     [st,flag,relres,iter,resvec]=gmres('AtimesX',Fine_bt, 960, tol, 1200,'Pre_D2_ichol',[], ini_sol);
%  
%     tEndD2 = cputime - tStart
%     flag
%     iter 
%     fprintf(fileID,'%5d %5d %6d\n', flag, iter);
% 
% %     Pre_D2=blkdiag(Fine_prob.Au, Fine_prob.SS_ksi, Fine_prob.SS_p);    %correct, (2.4) 02/04
% %     [eigvs, D, flg]=eigs(inv(Pre_D2)*Fine_prob.At, 654);
% %     eig_vals=diag(D)
% %     flg
% %     fprintf(fileID,'%5d \n', eig_vals);
% %     figure
% %     plot(eig_vals, 'o');    
% %     return;
%     
% 
%     fprintf(fileID, 'Block Diag with - sign in S_ksi \n');                % P_D3
%     tStart = cputime;
%     [st,flag,relres,iter,resvec]=gmres('AtimesX',Fine_bt, 960, tol, 1200,'Pre_D3_ichol',[], ini_sol);
% %    [st,flag,relres,iter ]=qmr('AtimesX',Fine_bt, tol, 1200,'Pre_D1_ichol', [] );
%     tEndD3 = cputime - tStart
%     flag
%     iter 
%     fprintf(fileID,'%5d %5d %6d\n', flag, iter);
% 
% %     Pre_D3=blkdiag(Fine_prob.Au, -Fine_prob.SS_ksi, -Fine_prob.SS_p);
% %     %eigs(inv(Pre_D3)*Fine_prob.At, 16)
% %     [eigvs, D, flg]=eigs(inv(Pre_D3)*Fine_prob.At, 654);
% %     eig_vals=diag(D)
% %     flg
% %     fprintf(fileID,'%5d \n', eig_vals);
% %     figure
% %     plot(eig_vals, '+');  
% %     return;
%     
%     fprintf(fileID, 'Block Diag with - sign in both S_ksi and S_p  \n');   % P_D4
%     tStart = cputime;
%     [st,flag,relres,iter,resvec]=gmres('AtimesX',Fine_bt, 960, tol, 1200,'Pre_D4_ichol',[], ini_sol);
% %    [st,flag,relres,iter ]=qmr('AtimesX',Fine_bt, tol, 1200,'Pre_D1_ichol', [] );
% 
%     tEndD4 = cputime - tStart
%     flag
%     iter 
%     fprintf(fileID,'%5d %5d %6d\n', flag, iter);
% %     Pre_D4=blkdiag(Fine_prob.Au, -Fine_prob.SS_ksi, Fine_prob.SS_p);
% %     [eigvs, D, flg]=eigs(inv(Pre_D4)*Fine_prob.At, 654);
% %     eig_vals=diag(D)
% %     flg
% %     fprintf(fileID,'%5d \n', eig_vals);
% %     figure
% %     plot(eig_vals, '*');
%     return;
    
    % Block Triangular preconditioners
    fprintf(fileID, 'Block triangular preconditioners, P_T1 \n');    
    % added by MCai, for block triangular case P1
    disp('P_T1');
    tStart = cputime;
    [st,flag,relres,iter,resvec]=gmres('AtimesX',Fine_bt, 960, tol, 1200,'Pre_T1_ichol',[], ini_sol);
%   [st,flag,relres,iter ]=qmr('AtimesX',Fine_bt, tol, 1200,'Pre_D1_ichol', [] );

	tEndT1 = cputime - tStart
    iter
    fprintf(fileID,'%5d %5d %6d\n', flag, iter);
   
%     Pre_T_P1=[Fine_prob.Au sparse(size(A_u, 1),size(SS_ksi, 2)) sparse(size(A_u, 1),size(SS_p, 2));
%               B_uksi  -Fine_prob.SS_ksi sparse(size(SS_ksi, 1),size(SS_p, 2));
%               sparse(size(SS_p, 1), size(A_u,2))  B_pksi    Fine_prob.SS_p];  % correct, MCai 02/04 (2.1)
%     [eigvs, D, flg]=eigs(inv(Pre_T_P1)*Fine_prob.At, 654);
%     eig_vals=diag(D)
%     flg
%     fprintf(fileID,'%5d \n', eig_vals);
%     figure
%     plot(eig_vals, '^r');   
    
    %MCai, for comparison put a negative sign
    fprintf(fileID,'P_T2');
    tStart = cputime;
    [st,flag,relres,iter,resvec]=gmres('AtimesX',Fine_bt, 960, tol, 1200,'Pre_T2_ichol',[], ini_sol);
%        [st,flag,relres,iter ]=qmr('AtimesX',Fine_bt, tol, 1200,'Pre_D1_ichol',[] );

	tEndT2 = cputime - tStart
    iter
    fprintf(fileID,'%5d %5d %6d\n', flag, iter);
%     Pre_T_M=blkdiag(Fine_prob.Au, -Fine_prob.SS_ksi, -Fine_prob.SS_p);
%     Pre_T_P2=[Fine_prob.Au sparse(size(A_u, 1),size(SS_ksi, 2)) sparse(size(A_u, 1),size(SS_p, 2));
%               B_uksi    Fine_prob.SS_ksi   sparse(size(SS_ksi, 1),size(SS_p, 2));
%               sparse(size(SS_p, 1), size(A_u,2)) -B_pksi     Fine_prob.SS_p];  % correct, MCai (2.3)
%     [eigvs, D, flg]=eigs(inv(Pre_T_P2)*Fine_prob.At, 654);
%     eig_vals=diag(D);
%     flg
%     fprintf(fileID,'%5d \n', eig_vals);
%     figure
%     plot(eig_vals, 'or');   
%     return;
    
    fprintf(fileID, 'P_T3');
    tStart = cputime;
    [st,flag,relres,iter,resvec]=gmres('AtimesX',Fine_bt, 960, tol, 1200,'Pre_T3_ichol',[], ini_sol);
 %     [st,flag,relres,iter ]=qmr('AtimesX',Fine_bt, tol, 1200,'Pre_D1_ichol', [] );

	tEndT3 = cputime - tStart
    iter
    fprintf(fileID,'%5d %5d %6d\n', flag, iter);

%     Pre_T_M=blkdiag(Fine_prob.Au, -Fine_prob.SS_ksi, -Fine_prob.SS_p);
%     Pre_T_P3=[Fine_prob.Au sparse(size(A_u, 1),size(SS_ksi, 2)) sparse(size(A_u, 1),size(SS_p, 2));
%               B_uksi  Fine_prob.SS_ksi sparse(size(SS_ksi, 1),size(SS_p, 2));
%               sparse(size(SS_p, 1), size(A_u,2))  -B_pksi    -Fine_prob.SS_p]; % correct, MCai ?2.3?SS_p ~= S_3
%     [eigvs, D, flg]=eigs(inv(Pre_T_P3)*Fine_prob.At, 654);
%      eig_vals=diag(D);
%     flg
%     fprintf(fileID,'%5d \n', eig_vals);
%     figure
%     plot(eig_vals, '+r');   
%     return;

    fprintf(fileID,'P_T4');
    tStart = cputime;
    [st,flag,relres,iter,resvec]=gmres('AtimesX',Fine_bt, 960, tol, 1200,'Pre_T4_ichol',[], ini_sol);
 %   [st,flag,relres,iter ]=qmr('AtimesX',Fine_bt, tol, 1200,'Pre_D1_ichol',[] );

	tEndT4 = cputime - tStart
    iter
     fprintf(fileID,'%5d %5d %6d\n', flag, iter);

%   %  Pre_T_M=blkdiag(Fine_prob.Au, -Fine_prob.SS_ksi, -Fine_prob.SS_p);
%     Pre_T_P4=[Fine_prob.Au sparse(size(A_u, 1),size(SS_ksi, 2)) sparse(size(A_u, 1),size(SS_p, 2));
%               B_uksi  -Fine_prob.SS_ksi sparse(size(SS_ksi, 1),size(SS_p, 2));
%               sparse(size(SS_p, 1), size(A_u,2))  B_pksi     -Fine_prob.SS_p];   %correct  ?2.3?
%     [eigvs, D, flg]=eigs(inv(Pre_T_P4)*Fine_prob.At, 654);
%      eig_vals=diag(D)
%     flg
%     fprintf(fileID,'%5d \n', eig_vals);
%     figure
%     plot(eig_vals, '*r');   
%     return;

%         relres_count(i_loop)=relres;
%         iter_count1(i_loop)=iter(1,2);
%         iter_count2(i_loop)=iter(1,1);
%         st=Fine_At\Fine_bt;
%         uvpnodes=[1:2*number_of_FE_nodes_u+2*number_of_FE_nodes_kp];
%         uvp_unkown=setdiff(uvpnodes,bound_uvp);
%         xst(uvp_unkown)=st;
%         ksi(:,1)=xst(1+2*number_of_FE_nodes_u:2*number_of_FE_nodes_u+number_of_FE_nodes_kp);
%         ph(:,1)=xst(1+2*number_of_FE_nodes_u+number_of_FE_nodes_kp:end);
    fprintf(fileID,'   The timestep is %d \n',timestep);
    fclose(fileID);
end


fprintf(1,'   Fine Biot done \n');    
    
%     uh1=xst(1:number_of_FE_nodes_u); 
%     [e_ul2,e_uh1]=errors_SU_new(element_area,T_basis_u,M_basis_u, uh1, number_of_elements, 6, 6,current_time); 
%     uh2=xst(1+number_of_FE_nodes_u:2*number_of_FE_nodes_u);
%     [e_vl2,e_vh1]=errors_SV_new(element_area,T_basis_u,M_basis_u, uh2, number_of_elements, 6, 6,current_time);
%     ksi=xst(1+2*number_of_FE_nodes_u:end);
%     [e_ksil2,e_ksih1]=errors_Sksi_new(element_area,T_partition,M_partition, ksi, number_of_elements, 3, 6,current_time);
%     [e_pl2,e_ph1]=errors_SP_new(element_area,T_partition,M_partition, ph, number_of_elements, 3, 6,current_time);
%      fprintf(1, '# FineNX =  :');       fprintf(1, '  %d \n', FineNx);   % which method we are using?
%     fprintf(1, '# e_ul2 =  :');        fprintf(1, '  %d \n', e_ul2);
%     fprintf(1, '# e_uh1 =  :');        fprintf(1, '  %d \n', e_uh1);
%    fprintf(1, '# e_vl2 =  :');        fprintf(1, '  %d \n', e_vl2);
%     fprintf(1, '# e_vh1 =  :');        fprintf(1, '  %d \n', e_vh1);
%     fprintf(1, '# e_ksil2 =  :');      fprintf(1, '  %d \n', e_ksil2);
%     fprintf(1, '# e_ksih1 =  :');      fprintf(1, '  %d \n', e_ksih1);
%     fprintf(1, '# e_pl2 =  :');        fprintf(1, '  %d \n', e_pl2);
%     fprintf(1, '# e_ph1 =  :');        fprintf(1, '  %d \n', e_ph1);
% 
% for i=1:4
%     rateuL2(i) = log(e_ul2(i)/e_ul2(i+1))/log(2.0);
% 	rateuH1(i) = log(e_uh1(i)/e_uh1(i+1))/log(2.0);
%     ratevL2(i) = log(e_vl2(i)/e_vl2(i+1))/log(2.0);
% 	ratevH1(i) = log(e_vh1(i)/e_vh1(i+1))/log(2.0);
% 	rateksiL2(i)=log(e_ksil2(i)/e_ksil2(i+1))/log(2.0);
% 	rateksiH1(i)=log(e_ksih1(i)/e_ksih1(i+1))/log(2.0);
% 	ratepL2(i)=log(e_pl2(i)/e_pl2(i+1))/log(2.0);
% 	ratepH1(i)=log(e_ph1(i)/e_ph1(i+1))/log(2.0);
%     fprintf(fid, '# rateuL2= :');      fprintf(fid, '  %d \n', rateuL2(i));
%      fprintf(fid, '# rateuH1= :');      fprintf(fid, '  %d \n', rateuH1(i));
%      fprintf(fid, '# rateksiL2= :');      fprintf(fid, '  %d \n', rateksiL2(i));
%      fprintf(fid, '# rateksiH1= :');      fprintf(fid, '  %d \n', rateksiH1(i));
%      fprintf(fid, '# ratepL2= :');      fprintf(fid, '  %d \n', ratepL2(i));
%      fprintf(fid, '# ratepH1= :');      fprintf(fid, '  %d \n', ratepH1(i));
% end
% 
% iter_count1
% relres_count
% iter_count2
% 
%       
% fprintf(fid, '# iter_count1= :');      fprintf(fid, '  %d \n', iter_count1);  % which method we are using?
% fprintf(fid, '# relres_count= :');      fprintf(fid, '  %d \n', relres_count);
% fprintf(fid, '# iter_count2= :');      fprintf(fid, '  %d \n', iter_count2);
% fprintf(fid,'KK=%d, alpha=%d, lambda=%d',[KK,alpha,lambda]);
% fprintf(fid,'########################################################## \n \n \n');
