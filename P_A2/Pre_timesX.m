function P_invX=Pre_timesX(rhs) 
% Actually, this part we only need to make the mapping property correct

global  Fine_prob 

P_invX=Fine_prob.Pt\rhs;

% P_invX=zeros(size(rhs));   % take 0 as initial solution 
% rhs_H=SubDom_prob(Nsubdomains+1).R*rhs;
% P_invX_C=SubDom_prob(Nsubdomains+1).At\rhs_H;
% 
% for i=1:Nsubdomains
%     rhs_i=SubDom_prob(i).R*rhs;
% 
%     u_i=SubDom_prob(i).At\rhs_i;
%     
%     P_invX=P_invX+(SubDom_prob(i).R)'*u_i;
% end
% P_invX=P_invX+(SubDom_prob(Nsubdomains+1).R)'*P_invX_C;
