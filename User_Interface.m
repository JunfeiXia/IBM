%% user interface version 2.0 March 15, 2020
%%%%%%%%%%%%%%%%%  input %%%%%%%%%%%%%%%%%
close all

year=100;    % year input

excel_output_index=0 % 1 means you want the excel output, it is everything 

%%% enter the figure number you want to plot, you can enter multiple numbers, like [1,2,10,11]...
figure_index=[1,2,3,4,5,6,7,8,9,10,11,12,13,14];
% Figure 1:  HEATMAP - Invasive leaf accumulation heat map
% Figure 2:  HEATMAP - Native leaf accumulation heat map
% Figure 3:  LINE CART - Invasive total trees trending
% Figure 4:  LINE CART - Native total trees trending
% Figure 5:  SCATTER - Total tree location
% Figure 6:  HISTOGRAM - Native trees d.b.h
% Figure 7:  HISTOGRAM - Invasive trees d.b.h
% Figure 8:  BAR - Number of new trees each year
% Figure 9:  BAR STACKED - Number of adult trees each year
% Figure 10:  BAR - Number of adult trees each year
% Figure 11: HEATMAP - Age of invasive trees
% Figure 12: HEATMAP - Age of native trees
% Figure 13: POLAR-HISTOGRAM # of Invasive trees and distance
% Figure 14: POLAR-HISTOGRAM # of Native trees and distance

%%%%%%%%%% clear %%%%%%%%%%%%%%%%
clear Tree_information_I
%%%%%%%%%%% format tranformation %%%%%%%%%%
% ID Type Age dbh xdot ydot
cnames={'ID','Type','Age','dbh','xdot','ydot'};
Tree_information_I(:,1)=cell2mat(Tree_information(year,1)); % ID
Tree_information_I(:,2)=cell2mat(Tree_information(year,2)); % Type
Tree_information_I(:,3)=cell2mat(Tree_information(year,3)); % Age
Tree_information_I(:,4)=cell2mat(Tree_information(year,4)); % dbh
Tree_information_I(:,5)=cell2mat(Tree_information(year,5)); % xdot
Tree_information_I(:,6)=cell2mat(Tree_information(year,6)); % ydot
% Tree_information_I(:,7)=cell2mat(Tree_information(year,7)); % Height
% Tree_information_I(:,8)=cell2mat(Tree_information(year,8)); % State
% Tree_information_I(:,9)=cell2mat(Tree_information(year,9)); % Number of years after dead


Tree_information_I_Invasive=Tree_information_I(Tree_information_I(:,2)==1,:);
Tree_information_I_Local=Tree_information_I(Tree_information_I(:,2)==2,:);

LF_annual_accumulation_invasive_I=cell2mat(LF_annual_accumulation_invasive(year));
LF_annual_accumulation_local_I=cell2mat(LF_annual_accumulation_local(year));

% data=dataset({Tree_information_I,'ID','Type','Age','dbh','xdot','ydot'});
clear IN_of_trees LN_of_trees
for i = 1:year
    temp = cell2mat(Tree_information(i,2));
    IN_of_trees(i)=length(find(temp==1));
    LN_of_trees(i)=length(find(temp==2));
end
%%%%%%% index setting %%%%%%%%%%
a_index=0;b_index=0;c_index=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%% display %%%%%%%%%%%%%
for display_i=1:length(figure_index)
    display_index=figure_index(display_i);
    switch display_index
        case 1
            figure(1) % invasive leaf accumulation heat map
            x_linespace=1:100;
            y_linespace=100:-1:1;
            pcolor(x_linespace,y_linespace,LF_annual_accumulation_invasive_I);
            colorbar
            title(['Invasive leaf accumulation for the ',num2str(year),'th year']);
            ylabel('grid 1:100');
            xlabel('grid 1:100');
        case 2
            figure(2)   % local leaf accumulation heat map
            % subplot(122)
            pcolor(x_linespace,y_linespace,LF_annual_accumulation_local_I);
            colorbar
            title(['Native leaf accumulation for the ',num2str(year),'th year']);
            ylabel('grid 1:100');
            xlabel('grid 1:100');

% I suggest to see the data in excel
% figure
% uitable('Data',Tree_information_I,'ColumnName',cnames);
% clear x_linespace y_linespace
        case 3
% scatter plot
            figure(3) % Invasive total trees trending
            x=[1:year]';
            % XI=[x,IN_of_trees'];XL=[x,LN_of_trees'];
            scatter(x,IN_of_trees);
            xlim([0,year]);%ylim([0,4000]);
            title(['invasive total trees trending ',num2str(year),'th year']);
            ylabel('tree total');
            xlabel('Nyear');
        case 4
            figure(4) % Native total trees trending
            scatter(x,LN_of_trees);
            xlim([0,year]);%ylim([0,1000]);
            title(['Native total trees trending ',num2str(year),'th year']);
            ylabel('tree total');
            xlabel('Nyear');
        case 5
            figure(5)  % Total tree location
            I=find(Tree_information_I(:,2)==1);
            scatter(Tree_information_I(I,5),Tree_information_I(I,6),'.','r'); % invasive tree location
            hold on;
            I=find(Tree_information_I(:,2)==2);
            scatter(Tree_information_I(I,5),Tree_information_I(I,6),'.','b');
            circle(50,50,5)
%             circle(50,50,10)
%             circle(50,50,20)
%             circle(50,50,30)
            title(['Total trees location ',num2str(year),'th year']);
            ylim([0,100]);ylabel('grid 1:100');
            xlim([0,100]);xlabel('grid 1:100');
        case 6
            figure(6) % Invasive trees d.b.h
            histogram(Tree_information_I_Invasive(:,4),[0:1:ceil(max(Tree_information_I_Invasive(:,4)))]);
            xlim([0,ceil(max(Tree_information_I(:,4)))]);
            title(['Invasive trees d.b.h at ',num2str(year),' year'])
            xlabel(['d.b.h (cm)']); ylabel('Number of trees');         
        case 7
            figure(7) % Native trees d.b.h
            histogram(Tree_information_I_Local(:,4),[0:1:ceil(max(Tree_information_I_Local(:,4)))]);
            xlim([0,ceil(max(Tree_information_I(:,4)))]);
            title(['Native trees d.b.h at ',num2str(year),' year'])
            xlabel(['d.b.h (cm)']); ylabel('Number of trees');
        case {8,9,10}
            if a_index==0
                figure(8) % Number of new trees each year
                clear n_NewTree_invasive n_NewTree_local
                for t=1:year
                    clear temp_Tree_information_I
                    temp_Tree_information_I(:,1)=cell2mat(Tree_information(t,2));
                    temp_Tree_information_I(:,2)=cell2mat(Tree_information(t,4));
                %     n_NewTree(t)=length(find(temp_Tree_information_I(:,2)==1.37));
                    temp_Tree_information_I_invasive = temp_Tree_information_I(find(temp_Tree_information_I(:,1)==1),2);
                    temp_Tree_information_I_local = temp_Tree_information_I(find(temp_Tree_information_I(:,1)==2),2);
                    n_trees_invasive(t)=length(find(temp_Tree_information_I_invasive>1.37));
                    n_trees_local(t)=length(find(temp_Tree_information_I_local>1.37));
                %     critical_ratio(t)=n_trees_invasive(t)/n_trees_local(t)
                    n_NewTree_invasive(t)=length(find(temp_Tree_information_I_invasive==1.37));
                    n_NewTree_local(t)=length(find(temp_Tree_information_I_local==1.37));
                end
                clear temp_Tree_information_I
                for i=1:(year/10)
                    Tenyear_n_NewTree_invasive(i)=sum(n_NewTree_invasive(i:10*i));
                    Tenyear_n_NewTree_local(i)=sum(n_NewTree_local(i:10*i));
                end
                n_NewTree=[n_NewTree_local',n_NewTree_invasive'];

                Tenyear_n_NewTree=[Tenyear_n_NewTree_local',Tenyear_n_NewTree_invasive'];
                bar(n_NewTree,'stacked');
                legend('native','invasive')
                title(['Number of saplings each year up to ',num2str(year),' year'])
                xlabel(['time in years']); ylabel('Number of trees');

                figure (9) % Number of adult trees each year

                n_trees=[n_trees_local',n_trees_invasive'];
                bar(n_trees,'stacked');
                legend('native','invasive')
                title(['Number of adult trees each year up to ',num2str(year),' year'])
                xlabel(['time in years']); ylabel('Number of trees');

                figure (10) % Number of adult trees each year
                bar(n_trees)
                legend('native','invasive')
                title(['Number of adult trees each year up to ',num2str(year),' year'])
                xlabel(['time in years']); ylabel('Number of trees');
                a_index=a_index+1;
            end

%%%%%%%%%%%%%%%% heat map age %%%%%%%%%%%%%%%%%%%%%%%%%
        case {11,12}
            if b_index==0
                x=0:100;
                y=100:-1:0;
                for i=1:100
                    for j=1:100
                        Inside=find(Tree_information_I(:,5)>x(i) & Tree_information_I(:,5)<x(i+1)...
                            & Tree_information_I(:,6)<y(j) & Tree_information_I(:,6)>y(j+1));
                        Inside_invasive=find(Tree_information_I(:,5)>x(i) & Tree_information_I(:,5)<x(i+1)...
                            & Tree_information_I(:,6)<y(j) & Tree_information_I(:,6)>y(j+1) & Tree_information_I(:,2)==1);
                        Inside_local=find(Tree_information_I(:,5)>x(i) & Tree_information_I(:,5)<x(i+1)...
                            & Tree_information_I(:,6)<y(j) & Tree_information_I(:,6)>y(j+1) & Tree_information_I(:,2)==2);

                        age_sum(j,i)=sum(Tree_information_I(Inside,3));
                        age_sum_invasive(j,i)=sum(Tree_information_I(Inside_invasive,3));
                        age_sum_local(j,i)=sum(Tree_information_I(Inside_local,3));
                        age_mean(j,i)=mean(Tree_information_I(Inside,3));
                    end
                end
                figure(11) % age of invasive
                pcolor(x_linespace,y_linespace,age_sum_invasive);
                colorbar;
                title(['age of invasive tree for ',num2str(year),' year']);
                ylabel('grid 1:100');
                xlabel('grid 1:100');

                figure(12) % age of local tree
                pcolor(x_linespace,y_linespace,age_sum_local);
                colorbar;
                title(['age of local tree for ',num2str(year),' year']);
                ylabel('grid 1:100');
                xlabel('grid 1:100');
                b_index=b_index+1;
            end
        case {13,14}
            if c_index==0
%%%%%%%% polar plot, density and distance %%%%%%%%%%%%%%%
                figure (13) % invasive
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
                polarhistogram(r_inside_circle_invasive*k,25);
                ax=gca;
                ax.ThetaTick = [0:36:360];
                ax.ThetaTickLabel={'0';'5';'10';'15';'20';'25';'30';'35';'40';'45';'50'};
                ax.RMinorTick = 'on';
                ax.ThetaMinorTick = 'on';
                title(['Invasive ', num2str(year), 'th year']);

                figure (14) % native
                polarhistogram(r_inside_circle_local*k,25);
                ax=gca;
                ax.ThetaTick = [0:36:360];
                ax.ThetaTickLabel={'0';'5';'10';'15';'20';'25';'30';'35';'40';'45';'50'};
                ax.RMinorTick = 'on';
                ax.ThetaMinorTick = 'on';
                title(['Native ', num2str(year), 'th year']);
                c_index=c_index+1;
            end
        otherwise
            '404 NOT FOUND'
    end
end
%% export to excel 
if excel_output_index==1
    excel_name=['Tree information for the ', num2str(year),'th year'];
    xlswrite([excel_name '.xlsx'],cnames,'sheet1','A1');
    xlswrite([excel_name '.xlsx'],Tree_information_I,'sheet1','A2');
end
%%%%%% clean the Workspace %%%%%%%%%%
clear excel_name a_index b_index c_index display_i display_index figure_index ax r_inside_circle_invasive r_inside_circle_local k 
clear temp_Tree_information_I temp_xdot_invasive temp_xdot_local temp_ydot_invasive temp_ydot_local temp_Tree_information_I_invasive temp_Tree_information_I_local
