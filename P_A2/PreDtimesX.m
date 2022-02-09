function P_invX=PreDtimesX(rhs) 
% Actually, this part we only need to make the mapping property correct

global  Nsubdomains SubDom_prob Coarse_prob level

if level==2

    P_invX=zeros(size(rhs));   % take 0 as initial solution
    rhs_H=SubDom_prob(Nsubdomains+1).R*rhs;
    P_invX_C=Coarse_prob.At\rhs_H;
    
    for i=1:Nsubdomains
        rhs_i=SubDom_prob(i).R*rhs;
        
        u_i=SubDom_prob(i).Pt\rhs_i;
        
        P_invX=P_invX+(SubDom_prob(i).R)'*u_i;
    end
    P_invX=P_invX+(SubDom_prob(Nsubdomains+1).R)'*P_invX_C;
else
    P_invX=zeros(size(rhs));   % take 0 as initial solution
    %rhs_H=SubDom_prob(Nsubdomains+1).R*rhs;
    %P_invX_C=Coarse_prob.At\rhs_H;
    
    for i=1:Nsubdomains
        rhs_i=SubDom_prob(i).R*rhs;
        
        u_i=SubDom_prob(i).Pt\rhs_i;
        
        P_invX=P_invX+(SubDom_prob(i).R)'*u_i;
    end
    %P_invX=P_invX+(SubDom_prob(Nsubdomains+1).R)'*P_invX_C;
end
    
    