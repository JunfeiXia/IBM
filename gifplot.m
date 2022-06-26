%%% visual gif %%%%
endyear=50;

export_index=1; % If you want to save the figures yealy, put 1 here. If you only want to see the gif, put 0 here.

delay_time=0.7; % the time between each figure

% put the location where you want to save the figures below
cd 'C:\Program Files\MATLAB\R2017a\bin\yuanming lu\2.28.2020\polar density and distance'


h=figure(1);
filename = 'Polarplot relationship between number of trees and distance.gif';

for year=1:endyear
    clear Tree_information_I
    %%%%%%%%%%% format tranformation section %%%%%%%%%%
    % ID Type Age dbh xdot ydot
    cnames={'ID','Type','Age','dbh','xdot','ydot'};
    Tree_information_I(:,1)=cell2mat(Tree_information(year,1));
    Tree_information_I(:,2)=cell2mat(Tree_information(year,2));
    Tree_information_I(:,3)=cell2mat(Tree_information(year,3));
    Tree_information_I(:,4)=cell2mat(Tree_information(year,4));
    Tree_information_I(:,5)=cell2mat(Tree_information(year,5));
    Tree_information_I(:,6)=cell2mat(Tree_information(year,6));
    % Tree_information_I(:,7)=cell2mat(Tree_information(year,7));
    LF_annual_accumulation_invasive_I=cell2mat(LF_annual_accumulation_invasive(year));
    LF_annual_accumulation_local_I=cell2mat(LF_annual_accumulation_local(year));
    
    clear IN_of_trees LN_of_trees
    for i = 1:year
        temp = cell2mat(Tree_information(i,2));
        IN_of_trees(i)=length(find(temp==1));
        LN_of_trees(i)=length(find(temp==2));
    end
    
    %%% copy the figure plot you want to see yearly from user interface and
    %%% paste it here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% plot section %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
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

k=2*pi/50;
r_inside_circle_invasive=rho_invasive(rho_invasive<=50);
r_inside_circle_local=rho_local(rho_local<=50);
polarhistogram(r_inside_circle_invasive*k,25);
ax=gca;
ax.ThetaTick = [0:36:360];
ax.ThetaTickLabel={'0';'5';'10';'15';'20';'25';'30';'35';'40';'45';'50'};
ax.RMinorTick = 'on';
ax.ThetaMinorTick = 'on';
title(['Invasive ', num2str(year), 'th year']);
%     pause (delay_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% gif export section %%%%%%%%%%%%%%%%%%%%%%%%%%
    drawnow
    % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if year == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delay_time); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delay_time); 
      end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% figure export section %%%%%%%%%%%%% 
    %%% save the figures one by one
    if export_index==1
        saveas(figure(1),['Invasive ',num2str(year), 'th year.png']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end