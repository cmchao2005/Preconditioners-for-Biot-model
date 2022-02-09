function Generate_Restriction(overlap_size, CoarseNx, CoarseNy, FineNx, FineNy,prob_lox, prob_loy, prob_hix, prob_hiy)
global Nsubdomains SubDom_prob level Fine_prob Coarse_prob enforce_local_0Diri detau
      
 
%use_Stokes_coarse=false;
SubDom_prob={};    % from 1 to Nsubdomains. 
% Find out the elements in Omega_i:  First figure out the boundary of each subdomains
i_sub=0;
for j=1:CoarseNy
  for i=1:CoarseNx
    i_sub=i_sub+1;
    if ((i-1)*Coarse_prob.hx-overlap_size+prob_lox > prob_lox-0.00001)
      SubDom_prob(i_sub).lbx=(i-1)*Coarse_prob.hx-overlap_size+prob_lox;
    else
      SubDom_prob(i_sub).lbx=(i-1)*Coarse_prob.hx+prob_lox;
    end
    if (i*Coarse_prob.hx+overlap_size+prob_lox < prob_hix+0.00001)
      SubDom_prob(i_sub).hbx=i*Coarse_prob.hx+overlap_size+prob_lox;
    else
      SubDom_prob(i_sub).hbx=i*Coarse_prob.hx+prob_lox;
    end 
    if ((j-1)*Coarse_prob.hy-overlap_size+prob_loy > prob_loy-0.00001)
      SubDom_prob(i_sub).lby=(j-1)*Coarse_prob.hy-overlap_size+prob_loy;
    else
      SubDom_prob(i_sub).lby=(j-1)*Coarse_prob.hy+prob_loy;
    end
    if (j*Coarse_prob.hy+overlap_size+prob_loy < prob_hiy+0.00001)
      SubDom_prob(i_sub).hby=j*Coarse_prob.hy+overlap_size+prob_loy;
    else
      SubDom_prob(i_sub).hby=j*Coarse_prob.hy+prob_loy;
    end
    
    % subdomain boundaries in non-overlapping sense 
    SubDom_prob(i_sub).lbx0=(i-1)*Coarse_prob.hx+prob_lox;
    SubDom_prob(i_sub).hbx0=i*Coarse_prob.hx+prob_lox;
    SubDom_prob(i_sub).lby0=(j-1)*Coarse_prob.hy+prob_loy;
    SubDom_prob(i_sub).hby0=j*Coarse_prob.hy+prob_loy;  
    
    % One layer shrinking
    if (SubDom_prob(i_sub).lbx > prob_lox)
      SubDom_prob(i_sub).lbx1=SubDom_prob(i_sub).lbx+Fine_prob.hx;
    else
      SubDom_prob(i_sub).lbx1=SubDom_prob(i_sub).lbx;
    end   
    if (SubDom_prob(i_sub).hbx < prob_hix)
      SubDom_prob(i_sub).hbx1=SubDom_prob(i_sub).hbx-Fine_prob.hx;
    else
      SubDom_prob(i_sub).hbx1=SubDom_prob(i_sub).hbx;
    end 
    if (SubDom_prob(i_sub).lby > prob_loy)
      SubDom_prob(i_sub).lby1=SubDom_prob(i_sub).lby+Fine_prob.hy;
    else
      SubDom_prob(i_sub).lby1=SubDom_prob(i_sub).lby;
    end
    if (SubDom_prob(i_sub).hby < prob_hiy)
      SubDom_prob(i_sub).hby1=SubDom_prob(i_sub).hby-Fine_prob.hy;
    else
      SubDom_prob(i_sub).hby1=SubDom_prob(i_sub).hby;
    end   
    
  end
end


for i_dom=1:Nsubdomains
  k=1; k1=1;
  for i_ele_c=1:2*FineNx*FineNy
    Center_point=Fine_prob.element_node_u(1:3, i_ele_c);
    Center_x=(Fine_prob.node_xy_u(1,Center_point(1))...
             +Fine_prob.node_xy_u(1,Center_point(2))...
             +Fine_prob.node_xy_u(1,Center_point(3)))/3;    
    Center_y=(Fine_prob.node_xy_u(2,Center_point(1))...
             +Fine_prob.node_xy_u(2,Center_point(2))...
             +Fine_prob.node_xy_u(2,Center_point(3)))/3;
    
    if ((Center_x >SubDom_prob(i_dom).lbx1) && (Center_y >SubDom_prob(i_dom).lby1) && ...
        (Center_x <SubDom_prob(i_dom).hbx1) && (Center_y <SubDom_prob(i_dom).hby1))
       SubDom_prob(i_dom).ELEs1(k1)=i_ele_c;  
       k1=k1+1; 
    end  
    
    if ((Center_x >SubDom_prob(i_dom).lbx) && (Center_y >SubDom_prob(i_dom).lby) && ...
        (Center_x <SubDom_prob(i_dom).hbx) && (Center_y <SubDom_prob(i_dom).hby))
      SubDom_prob(i_dom).elements(k)=i_ele_c;  
      k=k+1;
    end
  end
  
%   % Alternative way to determine the Boundary elements and the pressure nodes in outer layer
%   if (enforce_local_0Diri)
%     SubDom_prob(i_dom).Belements=setdiff(SubDom_prob(i_dom).elements, SubDom_prob(i_dom).ELEs1);
%   end   

  % determine the boundary nodes for each subdomain
  kkk=1;kkkp=1;   SubDom_prob(i_dom).bound={};
  for ele=1:size(SubDom_prob(i_dom).elements, 2)
    i_ele=SubDom_prob(i_dom).elements(ele);
    
    if abs(SubDom_prob(i_dom).lby-prob_loy)<0.001
        
        for i_node=1:6
          if ((abs(Fine_prob.node_xy_u(1,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).lbx) < 0.001*Fine_prob.hx) ||...
              (abs(Fine_prob.node_xy_u(1,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).hbx) < 0.001*Fine_prob.hx) ||...
              (abs(Fine_prob.node_xy_u(2,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).hby) < 0.001*Fine_prob.hx))
              ubound_node(i_dom, kkk)=Fine_prob.element_node_u(i_node,i_ele);
              kkk=kkk+1; 
          end

          if i_node<=3
              if ((abs(Fine_prob.node_xy_kp(1,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).lbx) < 0.001*Fine_prob.hx) ||...
                  (abs(Fine_prob.node_xy_kp(1,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).hbx) < 0.001*Fine_prob.hx) ||...
                  (abs(Fine_prob.node_xy_kp(2,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).hby) < 0.001*Fine_prob.hx))
                  pbound_node(i_dom, kkkp)=Fine_prob.element_node_kp(i_node,i_ele);
                  kkkp=kkkp+1; 
              end
          end
        end
        
    elseif abs(SubDom_prob(i_dom).hby-prob_hiy)<0.001
        
         for i_node=1:6
          if ((abs(Fine_prob.node_xy_u(1,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).lbx) < 0.001*Fine_prob.hx) ||...
              (abs(Fine_prob.node_xy_u(1,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).hbx) < 0.001*Fine_prob.hx) ||...
              (abs(Fine_prob.node_xy_u(2,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).lby) < 0.001*Fine_prob.hx)) 
              ubound_node(i_dom, kkk)=Fine_prob.element_node_u(i_node,i_ele);
              kkk=kkk+1; 
          end

          if i_node<=3
              if ((abs(Fine_prob.node_xy_kp(1,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).lbx) < 0.001*Fine_prob.hx) ||...
                  (abs(Fine_prob.node_xy_kp(1,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).hbx) < 0.001*Fine_prob.hx) ||...
                  (abs(Fine_prob.node_xy_kp(2,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).lby) < 0.001*Fine_prob.hx)) 
                  pbound_node(i_dom, kkkp)=Fine_prob.element_node_kp(i_node,i_ele);
                  kkkp=kkkp+1; 
              end
          end
         end
    else
         for i_node=1:6
          if ((abs(Fine_prob.node_xy_u(1,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).lbx) < 0.001*Fine_prob.hx) ||...
              (abs(Fine_prob.node_xy_u(1,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).hbx) < 0.001*Fine_prob.hx) ||...
              (abs(Fine_prob.node_xy_u(2,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).lby) < 0.001*Fine_prob.hx) ||...
              (abs(Fine_prob.node_xy_u(2,Fine_prob.element_node_u(i_node,i_ele))-SubDom_prob(i_dom).hby) < 0.001*Fine_prob.hx))
              ubound_node(i_dom, kkk)=Fine_prob.element_node_u(i_node,i_ele);
              kkk=kkk+1; 
          end

          if i_node<=3
              if ((abs(Fine_prob.node_xy_kp(1,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).lbx) < 0.001*Fine_prob.hx) ||...
                  (abs(Fine_prob.node_xy_kp(1,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).hbx) < 0.001*Fine_prob.hx) ||...
                  (abs(Fine_prob.node_xy_kp(2,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).lby) < 0.001*Fine_prob.hx) ||...
                  (abs(Fine_prob.node_xy_kp(2,Fine_prob.element_node_kp(i_node,i_ele))-SubDom_prob(i_dom).hby) < 0.001*Fine_prob.hx))
                  pbound_node(i_dom, kkkp)=Fine_prob.element_node_kp(i_node,i_ele);
                  kkkp=kkkp+1; 
              end
          end
         end
    end
  end
  
  utmp=ubound_node(i_dom, :);
  u_tmp_bound=utmp(utmp~=0);
  Bound_sort=sort(reshape(u_tmp_bound, [], 1));
  Bound_diff=diff(Bound_sort);
  Nodal_bound_u=[Bound_sort((Bound_diff~=0)); Bound_sort(end)];
  %Nodal_Dir_bound_u=setdiff(Nodal_bound_u,Fine_prob.Numan_bound_u);
  SubDom_prob(i_dom).bound_u=Nodal_bound_u; 
  
  ptmp=pbound_node(i_dom, :);
  p_tmp_bound=ptmp(ptmp~=0);
  Bound_sort_p=sort(reshape(p_tmp_bound, [], 1));
  Bound_diff_p=diff(Bound_sort_p);
  Nodal_bound_p=[Bound_sort_p((Bound_diff_p~=0)); Bound_sort_p(end)];
  %Nodal_Dir_bound_p=setdiff(Nodal_bound_p,Fine_prob.Numan_bound_p);
  SubDom_prob(i_dom).bound_p=Nodal_bound_p; 
  

  if (enforce_local_0Diri)
    
    temp_Nodal_ksi= Fine_prob.element_node_kp(:,SubDom_prob(i_dom).elements);
    Nodal_ksi=sort(reshape(temp_Nodal_ksi, [], 1));
    Node_ksi_diff=diff(Nodal_ksi);
    Nodal_Sksi=[Nodal_ksi((Node_ksi_diff~=0)); Nodal_ksi(end)]; 
    
    temp_Nodal_ksi_shink=Fine_prob.element_node_kp(:,SubDom_prob(i_dom).ELEs1);
    Nodal_ksi_shink=sort(reshape(temp_Nodal_ksi_shink, [], 1));
    Node_ksi_shink_diff=diff(Nodal_ksi_shink);
    Nodal_shink_Sksi=[Nodal_ksi_shink((Node_ksi_shink_diff~=0)); Nodal_ksi_shink(end)]; 
    ksiDiriNodes=setdiff(Nodal_Sksi,Nodal_shink_Sksi); % how about put more DOF, including derivatives.
    SubDom_prob(i_dom).ksiDiriNodes=ksiDiriNodes; 
  end
end

 

for i=1:Nsubdomains
  ele_i=SubDom_prob(i).elements;     % all the elements contained in the subdomain  
  Uindex=Fine_prob.element_node_u(1:6,ele_i);    % mark all nodes that in the subdomain ele_i
  
  % get the nodal index, delete those repeated nodal index
  Uindex_sort=sort(reshape(Uindex, [], 1));
  Udiff_index=diff(Uindex_sort);
  Nodal_Indx = [Uindex_sort((Udiff_index~=0)); Uindex_sort(end)];  %子单元所有节点

  Nu=size(Fine_prob.unknown_nodes_u,2);
  Nodal_Indx_copy=setdiff(Nodal_Indx, SubDom_prob(i).bound_u);  %  子单元去边界节点  
  Nnodes=size(Nodal_Indx_copy, 1);  %SubDom_prob(i).Nnodes=Nnodes;  % Only interior nodes and boundary nodes are kept
  SubDom_prob(i).Ru=sparse(1:Nnodes, Fine_prob.mv_node_u(Nodal_Indx_copy), ones(1, Nnodes), Nnodes,Nu);
  SubDom_prob(i).Nodal_Indx=Nodal_Indx;                       % keep this value  
  

  % restriction for p
  ksi_index=Fine_prob.element_node_kp(1:3,ele_i);
  ksi_index_sort=sort(reshape(ksi_index, [], 1));
  ksidiff_index=diff(ksi_index_sort);
  ksiNodal_Indx = [ksi_index_sort((ksidiff_index~=0)); ksi_index_sort(end)];
  
  Nksi=Fine_prob.number_of_node_kp;
  if (enforce_local_0Diri)
      ksiNodal_Indx_tot_mvc=setdiff(ksiNodal_Indx,SubDom_prob(i).ksiDiriNodes);
%       fileb=['arryp_2_2',mat2str(npart),mat2str(i),'.txt'];
%       dlmwrite(fileb,PNodal_Indx_tot_mvc,'delimiter','\t');
      ksiNnodes=size(ksiNodal_Indx_tot_mvc,1);
      Rksi_tmp=sparse(1:ksiNnodes, ksiNodal_Indx_tot_mvc,ones(ksiNnodes,1),ksiNnodes,Nksi);
    
  else
      ksiNnodes=size(ksiNodal_Indx,1);
      Rksi_tmp=sparse(1:ksiNnodes, ksiNodal_Indx,ones(ksiNnodes,1),ksiNnodes,Nksi);
  end
  SubDom_prob(i).Rksi=Rksi_tmp; 
  SubDom_prob(i).nksi_unknowns=ksiNnodes;
  
  Pindex=Fine_prob.element_node_kp(1:3,ele_i);    % mark all nodes that in the subdomain ele_i
  Pindex_sort=sort(reshape(Pindex, [], 1));
  Pdiff_index=diff(Pindex_sort);
  Nodal_Indx_P = [Pindex_sort((Pdiff_index~=0)); Pindex_sort(end)];  %子单元所有节点

  Np=size(Fine_prob.unknown_nodes_p,2);
  Nodal_Indx_copy_p=setdiff(Nodal_Indx_P, SubDom_prob(i).bound_p);  %  子单元去边界节点  
  Npnodes=size(Nodal_Indx_copy_p, 1);  %SubDom_prob(i).Nnodes=Nnodes;  % Only interior nodes and boundary nodes are kept
  SubDom_prob(i).Rp=sparse(1:Npnodes, Fine_prob.mv_node_p(Nodal_Indx_copy_p), ones(1, Npnodes), Npnodes,Np);
  

  R = [SubDom_prob(i).Ru     sparse(Nnodes, Nu)   sparse(Nnodes,Nksi)   sparse(Nnodes,Np);
        sparse(Nnodes, Nu)   SubDom_prob(i).Ru     sparse(Nnodes,Nksi)   sparse(Nnodes,Np);
        sparse(ksiNnodes, Nu)  sparse(ksiNnodes,Nu)   SubDom_prob(i).Rksi   sparse(ksiNnodes,Np);
        sparse(Npnodes, Nu)        sparse(Npnodes,Nu)         sparse(Npnodes,Nksi)    SubDom_prob(i).Rp]; 

  SubDom_prob(i).R=sparse(R);
  
%   u_1=ones(size(SubDom_prob(i).Rksi,1),1);
%   uuu=SubDom_prob(i).Rksi'*u_1;
%   BBBu=reshape(uuu',[FineNx+1,FineNx+1]);
%   BBBu'
%    [iR,jR,valuR] = find(SubDom_prob(i).R);
%    data_dump = [iR,jR,valuR];
%    fileu=['R_',mat2str(CoarseNx),'_',mat2str(detau),'_',mat2str(i),'_m.txt'];
%    fid = fopen(fileu,'w');
%    fprintf( fid,'%d %d %f\n', transpose(data_dump) );
%    fclose(fid);
  
    %   Ruv = [SubDom_prob(i).Ru sparse(Nnodes, Nu);
    %         sparse(Nnodes, Nu) SubDom_prob(i).Rv]; 
    %   SubDom_prob(i).Ruv=sparse(Ruv);
    %   
    %   SubDom_prob(i).At=SubDom_prob(i).R*Fine_prob.At*(SubDom_prob(i).R).';   % build up the local matrix
    SubDom_prob(i).Pt=SubDom_prob(i).R*Fine_prob.Pt*(SubDom_prob(i).R).';
    %   SubDom_prob(i).Auv=SubDom_prob(i).Ruv*Fine_prob.Auv*(SubDom_prob(i).Ruv');


end

if (level==2) 
  CNu_nodes=(2*CoarseNx+1)*(2*CoarseNy+1);
  CNkp_nodes=(CoarseNx+1)*(CoarseNy+1);     % equal to number of elements 
  Nu_nodes=(2*FineNx+1)*(2*FineNy+1);
  Nkp_nodes=(FineNx+1)*(FineNy+1);      
  RH_u=zeros(CNu_nodes, Nu_nodes);       % RH_u is from fine to coarse 
  RH_ksi=zeros(CNkp_nodes, Nkp_nodes);
  RH_p=zeros(CNkp_nodes, Nkp_nodes);
   
  for i=1:Nu_nodes
     % first get the coordinates of i_th node on fine mesh. 
     i_xcoor=Fine_prob.node_xy_u(1,i);
     i_ycoor=Fine_prob.node_xy_u(2,i); 

     for i_ele=1:size(Coarse_prob.element_node_u, 2)

        % get the 4 corner points of a coarse element
        X_vtx(1:3,1)=Coarse_prob.node_xy_u(1,Coarse_prob.element_node_u(1:3,i_ele));
        Y_vtx(1:3,1)=Coarse_prob.node_xy_u(2,Coarse_prob.element_node_u(1:3,i_ele));  
        
        Ae=[X_vtx(2)-X_vtx(1), X_vtx(3)-X_vtx(1);
            Y_vtx(2)-Y_vtx(1), Y_vtx(3)-Y_vtx(1)];
        be=[i_xcoor-X_vtx(1);i_ycoor-Y_vtx(1)];
        ue=Ae\be;


        Ae1=[X_vtx(3)-X_vtx(2), X_vtx(1)-X_vtx(2);
             Y_vtx(3)-Y_vtx(2), Y_vtx(1)-Y_vtx(2)];
        be1=[i_xcoor-X_vtx(2);i_ycoor-Y_vtx(2)];
        ue1=Ae1\be1;
        
        Ae2=[X_vtx(1)-X_vtx(3), X_vtx(2)-X_vtx(3);
             Y_vtx(1)-Y_vtx(3), Y_vtx(2)-Y_vtx(3)];
        be2=[i_xcoor-X_vtx(3);i_ycoor-Y_vtx(3)];
        ue2=Ae2\be2;
    
        if ( (ue(1) >=-0.00001) && (ue(1) <=1.00001) && (ue(2) >=-0.00001) && (ue(2) <=1.00001) &&...
             (ue1(1) >=-0.00001) && (ue1(1) <=1.00001) && (ue1(2) >=-0.00001) && (ue1(2) <=1.00001)&& ...
             (ue2(1) >=-0.00001) && (ue2(1) <=1.00001) && (ue2(2) >=-0.00001) && (ue2(2) <=1.00001))
            for j=1:6  
                in=Coarse_prob.element_node_u(j,i_ele );
                RH_u(in, i)=u_basis(j,ue(1),ue(2));
            end 
        end
        
     end
  end

 RH_u(:, Fine_prob.bound_u)=[];
 RH_u(Coarse_prob.bound_u, :)=[];
  
  for i=1:Nkp_nodes
        % first get the coordinates of i_th node on fine mesh. 
        i_xcoor_kp=Fine_prob.node_xy_kp(1,i);
        i_ycoor_kp=Fine_prob.node_xy_kp(2,i);
     
      for i_ele_kp=1:size(Coarse_prob.element_node_kp, 2)
          
        % get the 4 corner points of a coarse element
        Xkp_vtx(1:3,1)=Coarse_prob.node_xy_kp(1,Coarse_prob.element_node_kp(1:3,i_ele_kp));
        Ykp_vtx(1:3,1)=Coarse_prob.node_xy_kp(2,Coarse_prob.element_node_kp(1:3,i_ele_kp));  
      
        Ape=[Xkp_vtx(2)-Xkp_vtx(1), Xkp_vtx(3)-Xkp_vtx(1);
            Ykp_vtx(2)-Ykp_vtx(1), Ykp_vtx(3)-Ykp_vtx(1)];
        bpe=[i_xcoor_kp-Xkp_vtx(1);i_ycoor_kp-Ykp_vtx(1)];
        pe=Ape\bpe;


        Ape1=[Xkp_vtx(3)-Xkp_vtx(2), Xkp_vtx(1)-Xkp_vtx(2);
             Ykp_vtx(3)-Ykp_vtx(2), Ykp_vtx(1)-Ykp_vtx(2)];
        bpe1=[i_xcoor_kp-Xkp_vtx(2);i_ycoor_kp-Ykp_vtx(2)];
        pe1=Ape1\bpe1;
        
        Ape2=[Xkp_vtx(2)-Xkp_vtx(3), Xkp_vtx(1)-Xkp_vtx(3);
             Ykp_vtx(2)-Ykp_vtx(3), Ykp_vtx(1)-Ykp_vtx(3)];
        bpe2=[i_xcoor_kp-Xkp_vtx(3);i_ycoor_kp-Ykp_vtx(3)];
        pe2=Ape2\bpe2;

  
        if ( (pe(1) >=-0.00001) && (pe(1) <=1.00001) && (pe(2) >=-0.00001) && (pe(2) <=1.00001) &&...
             (pe1(1) >=-0.00001) && (pe1(1) <=1.00001) && (pe1(2) >=-0.00001) && (pe1(2) <=1.00001)&& ...
             (pe2(1) >=-0.00001) && (pe2(1) <=1.00001) && (pe2(2) >=-0.00001) && (pe2(2) <=1.00001))

            for j=1:3
                in=Coarse_prob.element_node_kp(j,i_ele_kp );
                RH_ksi(in, i)=p_basis(j,pe(1),pe(2));
            end 
         end  
     end   
  end

  nu_unknowns_coarse=CNu_nodes-size(Coarse_prob.bound_u,2);
  np_unknowns_coarse=CNkp_nodes-size(Coarse_prob.bound_p,2);
  nu_unknowns_fine=Fine_prob.number_of_unknown_u;
  np_unknowns_fine=Fine_prob.number_of_unknown_p;
  RH_p=RH_ksi;
  RH_p(:, Fine_prob.bound_p)=[]; RH_p( Coarse_prob.bound_p,:)=[];
  R0=[RH_u      sparse(nu_unknowns_coarse, nu_unknowns_fine)          sparse(nu_unknowns_coarse, Fine_prob.number_of_node_kp)   sparse(nu_unknowns_coarse,np_unknowns_fine);
      sparse(nu_unknowns_coarse, nu_unknowns_fine)      RH_u          sparse(nu_unknowns_coarse, Fine_prob.number_of_node_kp)   sparse(nu_unknowns_coarse,np_unknowns_fine);
      sparse(Coarse_prob.number_of_node_kp, nu_unknowns_fine) sparse(Coarse_prob.number_of_node_kp,  nu_unknowns_fine)    RH_ksi   sparse(Coarse_prob.number_of_node_kp,np_unknowns_fine);
      sparse(np_unknowns_coarse, nu_unknowns_fine)                            sparse(np_unknowns_coarse,  nu_unknowns_fine)    sparse(np_unknowns_coarse, Fine_prob.number_of_node_kp)     RH_p ];

  SubDom_prob(Nsubdomains+1).R=sparse(R0);
  
  SubDom_prob(Nsubdomains+1).Rksi=sparse(RH_ksi);
  SubDom_prob(Nsubdomains+1).Rp=sparse(RH_p);
  
%    [i0,j0,val0] = find( SubDom_prob(Nsubdomains+1).R);
%    data_dump0 = [i0,j0,val0];
%    fileR0=['R0_',mat2str(CoarseNx),'_m.txt'];
%    fidr0 = fopen(fileR0,'w');
%    fprintf( fidr0,'%d %d %f\n', transpose(data_dump0) );
%    fclose(fidr0);
                    
  
end


