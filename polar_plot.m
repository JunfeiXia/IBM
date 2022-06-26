%%%%%  polarplot  %%%%
close all

year=1;    % year input

%%%%%%%%%% clear %%%%%%%%%%%%%%%%
clear Tree_information_I

%%%%%%% read the tree information data %%%%%%
% ID Type Age dbh xdot ydot
cnames={'ID','Type','Age','dbh','xdot','ydot'};
Tree_information_I(:,1)=cell2mat(Tree_information(year,1));
Tree_information_I(:,2)=cell2mat(Tree_information(year,2));
Tree_information_I(:,3)=cell2mat(Tree_information(year,3));
Tree_information_I(:,4)=cell2mat(Tree_information(year,4));
Tree_information_I(:,5)=cell2mat(Tree_information(year,5));
Tree_information_I(:,6)=cell2mat(Tree_information(year,6));
% Tree_information_I(:,7)=cell2mat(Tree_information(year,7));

%%%%%%%% transform x-y coordinate to polar coordinate %%%%%
%%%% step 1, change the original point from (0,0) to (50,50)
Tree_information_I_Invasive=Tree_information_I(Tree_information_I(:,2)==1,:);
Tree_information_I_Local=Tree_information_I(Tree_information_I(:,2)==2,:);
temp_xdot_invasive=Tree_information_I_Invasive(:,5)-50;
temp_ydot_invasive=Tree_information_I_Invasive(:,6)-50;
temp_xdot_local=Tree_information_I_Local(:,5)-50;
temp_ydot_local=Tree_information_I_Local(:,6)-50;
%%%% step 2, transform x-y to ploar
[theta_invasive,rho_invasive] = cart2pol(temp_xdot_invasive,temp_ydot_invasive);
Location_invasive=[theta_invasive,rho_invasive];
[theta_local,rho_local] = cart2pol(temp_xdot_local,temp_ydot_local);
Location_local=[theta_local,rho_local] ;
%%%% step 3 plot %%% turn distance to degree, only consider the maximum circle inside the square 
% figure
% polarscatter(theta,rho,'.');
k=2*pi/50;
r_inside_circle_invasive=rho_invasive(rho_invasive<=50);
r_inside_circle_local=rho_local(rho_local<=50);
figure % invaisve
polarhistogram(r_inside_circle_invasive*k,25);
ax=gca;
ax.ThetaTick = [0:36:360];
ax.ThetaTickLabel={'0';'5';'10';'15';'20';'25';'30';'35';'40';'45';'50'};
ax.RMinorTick = 'on';
ax.ThetaMinorTick = 'on';
title(['Invasive ', num2str(year), 'th year']);
% Below is the formal title,use the simplified title first
% title(['Relationship between the number of invasive tree and the distance to the center for the ', num2str(year), 'th year']);
% polarhistogram(theta,20); % density & direction, you might need this figure




%%%%%% density %%%%%%% abandoned
% theta_linespace=[-pi:pi/2:pi];
% rho_linespace=[-50:1:50];
% theta_linecenter=[0.5:1:359.5];
% rho_linecenter=[-49.5:1:49.5];
% clear density_Z_invasive density_Z_local
% for theta_i=1:360
%     for r_j=1:100
%         density_Z_invasive(theta_i,r_j)=length(find(Location_invasive(:,1)>=theta_linespace(theta_i) & Location_invasive(:,1)<theta_linespace(theta_i+1) ...
%             & Location_invasive(:,2)>=rho_linespace(r_j) &  Location_invasive(:,2)<rho_linespace(r_j+1)));
%         density_Z_local(theta_i,r_j)=length(find(Location_local(:,1)>=theta_linespace(theta_i) & Location_local(:,1)<theta_linespace(theta_i+1) ...
%             & Location_local(:,2)>=rho_linespace(r_j) &  Location_local(:,2)<rho_linespace(r_j+1)));
%     end
% end
    
