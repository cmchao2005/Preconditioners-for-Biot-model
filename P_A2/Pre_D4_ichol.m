function P_invX=Pre_D4_ichol(rhs) 
% Actually, this part we only need to make the mapping property correct

global  Fine_prob chol_Au chol_SSksi chol_SSp


rhs1=rhs(1:size(Fine_prob.Au,1));
rhs2=rhs(1+size(Fine_prob.Au,1):size(Fine_prob.Au,1)+size(Fine_prob.SS_ksi,1));
rhs3=rhs(1+size(Fine_prob.Au,1)+size(Fine_prob.SS_ksi,1):end);

%disp('pre_d is called');

%chol_Au= ichol(Fine_prob.Au);
P_invX_u= (chol_Au.')\(chol_Au\rhs1);
%chol_SSksi= ichol(Fine_prob.SS_ksi);
P_invX_ksi=-(chol_SSksi.')\(chol_SSksi\rhs2);
%chol_SSp= ichol(-Fine_prob.SS_p);
P_invX_p=(chol_SSp.')\(chol_SSp\rhs3);
 
P_invX=[P_invX_u;P_invX_ksi;P_invX_p];
 