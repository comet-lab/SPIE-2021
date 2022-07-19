close all;

exp_left_simID = 'CoronalView2-left-nowrist-dq-0-91pts';
exp_right_simID = 'CoronalView2-right-nowrist-dq-0-91pts';
exp_tot_simID = 'CoronalView2-nowrist-dq-0-91pts';
left_simID = 'CoronalView2-left-nowrist-dq-0.06-10000pts';

reg_c = [1 0.5 0.45];
vis_c = [47,19,16]/255; % dark purple
vis_c = [89 248 255]/255
color_map = [reg_c; vis_c];

%%
figure
% subplot(121)
makeVisibilityFig(exp_left_simID,'plotVisible',true,'plotCones',false,'tissueName','L1-surface-op4-left','colorMap',color_map);
xlabel('X [mm]','FontSize',14,'FontName','CMU Serif');
ylabel('Y [mm]','FontSize',14,'FontName','CMU Serif');
zlabel('Z [mm]','FontSize',14,'FontName','CMU Serif');
title('','FontSize',16,'FontName','CMU Serif')
view(70,20)
curr_lims = zlim();
zlim([curr_lims(1),-25])

% subplot(122)
figure
makeVisibilityFig(exp_right_simID,'plotVisible',true,'plotCones',false, 'tissueName','L1-surface-op4-right','colorMap',color_map);
xlabel('X [mm]','FontSize',14,'FontName','CMU Serif');
ylabel('Y [mm]','FontSize',14,'FontName','CMU Serif');
zlabel('Z [mm]','FontSize',14,'FontName','CMU Serif');
title('','FontSize',16,'FontName','CMU Serif')
curr_lims = zlim();
zlim([curr_lims(1),-25])
view(-70,20)

%%
close all
figure
makeVisibilityFig(exp_tot_simID,'plotVisible',true,'plotCones',false,'tissueName','tissue_cropped','colorMap',color_map);
xlabel('X [mm]','FontSize',14,'FontName','CMU Serif');
ylabel('Y [mm]','FontSize',14,'FontName','CMU Serif');
zlabel('Z [mm]','FontSize',14,'FontName','CMU Serif');
title('','FontSize',16,'FontName','CMU Serif')
view(-50,20)
curr_lims = zlim();
title('')
% zlim([curr_lims(1),-25])
axis off
grid off

%%
% Right Vocal Fold 
figure
makeVisibilityFig(exp_tot_simID,'plotVisible',true,'plotCones',false,'tissueName','tissue_cropped','colorMap',color_map);
xlabel('X [mm]','FontSize',14,'FontName','CMU Serif');
ylabel('Y [mm]','FontSize',14,'FontName','CMU Serif');
zlabel('Z [mm]','FontSize',14,'FontName','CMU Serif');
title('','FontSize',16,'FontName','CMU Serif')
view(-50,20)
curr_lims = zlim();
zlim([-50,-30])
xlim([28,36])
ylim([-26,-21])


% Left Vocal Fold
figure
makeVisibilityFig(exp_tot_simID,'plotVisible',true,'plotCones',false,'tissueName','tissue_cropped');
xlabel('X [mm]','FontSize',14,'FontName','CMU Serif');
ylabel('Y [mm]','FontSize',14,'FontName','CMU Serif');
zlabel('Z [mm]','FontSize',14,'FontName','CMU Serif');
title('','FontSize',16,'FontName','CMU Serif')
view(50,20)
curr_lims = zlim();
zlim([-50,-30])
xlim([21,28])
ylim([-26,-21])