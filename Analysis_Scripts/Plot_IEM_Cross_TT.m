%{
Plot_IEM_Cross_TT
Author: Tom Bullock, UCSB Attention Lab
Date: 11.18.19

%}

clear 
close all

% which directories?
rDir = '/Users/tombullock/Documents/Psychology/WTF_EYE';
sourceDir = [rDir '/' 'IEM_Results_TT_Cross'];
destDir = [rDir '/' 'Data_Compiled'];

% which subjects?
subs = [5];

% compile slope data
for iSub=1:length(subs)
    
    sjNum=subs(iSub);
    
    % load data
    load([sourceDir '/' sprintf('sj%02d_IEM_Alpha_CrossTT.mat',sjNum)])
    allTF_real_total(iSub,:,:,:)=em.rSl.total; %subs x cross x tr x te
    clear em 
end


xNew=-500:2501/40:2000;
yNew=xNew;
%cLims = [-.0003 .0007];
%cLims = [0 .0007]

h=figure;

for i=1:8 %size(allTF_real_total,2)

    if i==2; tl='Eyes Closed'; m=13;
    elseif i==1; tl='Eyes Open'; m=14;
    elseif i==4; tl='Eyes Closed Masked'; m=15;
    elseif i==3; tl='Eyes Open Masked'; m=16;
    elseif i==6; tl='Train EC, Test EO'; m=1;
    elseif i==5; tl='Train EO, Test EC'; m=4;
    elseif i==8; tl='Train ECM, Test EOM'; m=9;
    elseif i==7; tl='Train EOM, Test ECM'; m=12;
    end
    
    subplot(2,4,i)
    imagesc(xNew,yNew,squeeze(mean(allTF_real_total(:,m,:,:),1))) %cLims
    ylabel('Training')
    xlabel('Testing')
    pbaspect([1,1,1])
    title(tl)
    cbar
    
    vline(0,'k--')
    hline(0,'k--')
    vline(250,'k--')
    hline(250,'k--')
    
    if i==3||i==4||i==7||i==8
    vline(300,'k--')
    hline(300,'k--')
    end
    %pause(1)
    
end

% isolate 8 important comparisons (see above plot)
%allTF_mat = (allTF_real_total(:,[13,14,15,16,1,4,9 12],:,:));

parsave([destDir '/' 'ALL_TT_WITHIN_BETWEEN.mat'],allTF_real_total,allTF_perm_total)

%% DISPLAY ALL 12 additional training/testing cross plots
h=figure;
for i=1:12
    
    if i==1; tl='TR-EC, TE-EO';
    elseif i==2; tl='TR-EC, TE-ECM';
    elseif i==3; tl='TR-EC, TE-EOM';
    elseif i==4; tl='TR-EO, TE-EC';
    elseif i==5; tl='TR-EO, TE-ECM';
    elseif i==6; tl='TR-EO, TE-EOM';
    elseif i==7; tl='TR-ECM, TE-EC';
    elseif i==8; tl='TR-ECM, TE-EO';
    elseif i==9; tl='TR-ECM, TE-EOM';
    elseif i==10; tl='TR-EOM, TE-EC';
    elseif i==11; tl='TR-EOM, TE-EO';
    elseif i==12; tl='TR-EOM, TE-ECM';
    end 

    subplot(4,3,i)
    imagesc(xNew,yNew,squeeze(mean(allTF_real_total(:,i,:,:),1))) %cLims
    ylabel('Training')
    xlabel('Testing')
    pbaspect([1,1,1])
    title(tl)
    %cbar
    vline(0,'k--')
    hline(0,'k--')
    vline(250,'k--')
    hline(250,'k--')
   
end