%{
Plot_IEM_Within
Author: Tom Bullock, UCSB Attention Lab
Date: 11.19.19

Plot TF surface plots.
%}

clear
close all

% load 
rDir = '/home/waldrop/Desktop/WTF_EYE';
sourceDir = [rDir '/' 'IEM_Results_TT_Within_Fixed'];
plotDir = '/home/waldrop/Desktop/WTF_EYE/Plots';

% subjects [TOM ADDED SUBJECT 7 BACK IN}!
subjects = [1:7,9:14,16:20,22:27,31];
%subjects = 31; % good strong response - use for eye-data?

% set up plot
thisSR=256;
xNewTick = [1, thisSR*.5, thisSR*1,  thisSR*1.5, thisSR*2, thisSR*2.5];
xNewTickLabel = [];%['-0.5';'  0 ';' 0.5';'  1 ';' 1.5';'  2 '];
xNewLim = [1 640];
yNewTick = [1 2 3 4 5 6 7 8];
zNewTick = [.2 .4 .6];
zNewTickLabel = [];%zNewTick;
yNewTickLabel = [];%['-180';'    ';'-90 ';'    ';' 0  ';'    ';' 90 ';'    '];
yNewLim = [1 8];
colorBar = [0 .8];
thisAxisFontsize = 20;

% % disable tick labels
% xNewTickLabel = [];
% yNewTickLabel = [];
% zNewTickLabel = [];

% compile data
for iSub=1:length(subjects)
    sjNum=subjects(iSub);
    load([sourceDir '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)],'em')
    
    all_real_total(iSub,1,:,:) = em.tf_total_cond1;
    all_real_total(iSub,2,:,:) = em.tf_total_cond2;
    all_real_total(iSub,3,:,:) = em.tf_total_cond3;
    all_real_total(iSub,4,:,:) = em.tf_total_cond4;

end

% plot data
f=1; % alpha, theta
for iCond=1:4
    
    h=figure('units','normalized','outerposition',[0 0 .65 1]);

    if      iCond==1; mc='Spatial'; ec='Fix';
    elseif  iCond==2; mc='Color';   ec='Fix';
    elseif  iCond==3; mc='Spatial'; ec='Move';
    elseif  iCond==4; mc='Color';   ec='Move';
    end
    
    %subplot(2,2,iCond);
    surf(squeeze(mean(all_real_total(:,iCond,:,:),1)),'linestyle','none','FaceAlpha',1); 

    %ylabel('Time(s)')
    %zlabel('Chan. Resp.(uV2)')
    title([mc ' / ' ec])
    view(70,30)
    zlim([0 .65])
    colormap('jet')
    caxis([.2,.6]);
    pbaspect([1 2 1])
    set(gca,'yTick',xNewTick,'yticklabel',xNewTickLabel,'xTick',yNewTick,...
        'xticklabel',yNewTickLabel,'ylim',xNewLim,'fontsize',thisAxisFontsize,...
        'LineWidth',4,'ztick',zNewTick,'zticklabel',zNewTickLabel);
    grid('off')
    
    saveas(h,[plotDir '/' 'IEM_Surf_Fixed' mc '_' ec '.eps'],'epsc')
    saveas(h,[plotDir '/' 'IEM_Surf_Fixed' mc '_' ec '.jpg'],'jpeg')
end

