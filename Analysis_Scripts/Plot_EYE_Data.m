clear
close all

% set dirs
sourceDir = 'D:\\WTF_EYE\\EYE\\Synchronized_EYE';

%subjects
sjNum=11;

% load data
load(sprintf([sourceDir '\\sj%02d_eye_beh_sync.mat'],sjNum))

% get example moving condition (cat session 1&2)
     %Location
trial_data1 = masterStruct(1,3).allTrialData;
trial_data2 = masterStruct(2,3).allTrialData;
trial_cat = [trial_data1,trial_data2];
     %Color
% trial_data1 = masterStruct(1,4).allTrialData;
% trial_data2 = masterStruct(2,4).allTrialData;
% trial_cat = [trial_data1,trial_data2];

%get rid of brokentrials
cnt=0;
for i=1:length(trial_cat)
    if trial_cat(i).brokenTrial
        cnt=cnt+1;
        brokenTrialVec(cnt)=i;
    end
end
trial_cat(brokenTrialVec) = [];

% cat session 1 & 2 eye data
     %Location
eye_dataX1 = masterStruct(1,3).eyeData.allXpos;
eye_dataX2 = masterStruct(2,3).eyeData.allXpos;
eye_catX = [eye_dataX1;eye_dataX2];
eye_dataY1 = masterStruct(1,3).eyeData.allYpos;
eye_dataY2 = masterStruct(2,3).eyeData.allYpos;
eye_catY = [eye_dataY1;eye_dataY2];
     %Color
% eye_dataX1 = masterStruct(1,4).eyeData.allXpos;
% eye_dataX2 = masterStruct(2,4).eyeData.allXpos;
% eye_catX = [eye_dataX1;eye_dataX2];
% eye_dataY1 = masterStruct(1,4).eyeData.allYpos;
% eye_dataY2 = masterStruct(2,4).eyeData.allYpos;
% eye_catY = [eye_dataY1;eye_dataY2];


% remove broken trials from eye_data/cat struct
%eye_data.allTimes(brokenTrialVec,:) = [];
eye_catX(brokenTrialVec,:) = [];
eye_catY(brokenTrialVec,:) = [];
%eye_data.allPA(brokenTrialVec,:) = [];
%eye_data.targCode(brokenTrialVec,:) = [];


for i=1:size(eye_catX,1)
    for j=1:size(eye_catX,2)
        if eye_catX(i,j)>2000
            eye_catX(i,j) = NaN;
        end
    end
end
for i=1:size(eye_catY,1)
    for j=1:size(eye_catY,2)
        if eye_catY(i,j)>2000
            eye_catY(i,j) = NaN;
        end
    end
end

%find eye_Loc trials
loc1 = []; loc1_count = 1;
loc2 = []; loc2_count = 1;
loc3 = []; loc3_count = 1;
loc4 = []; loc4_count = 1;
loc5 = []; loc5_count = 1;
loc6 = []; loc6_count = 1;
loc7 = []; loc7_count = 1;
loc8 = []; loc8_count = 1;

for iTrial=1:length(trial_cat)
    
    eye_Loc = trial_cat(iTrial).eyeLoc;
    if eye_Loc==1
        loc1(loc1_count) = iTrial;
        loc1_count = loc1_count + 1;
    elseif eye_Loc==2
        loc2(loc2_count) = iTrial;
        loc2_count = loc2_count + 1;
    elseif eye_Loc==3
        loc3(loc3_count) = iTrial;
        loc3_count = loc3_count + 1;
    elseif eye_Loc==4
        loc4(loc4_count) = iTrial;
        loc4_count = loc4_count + 1;
    elseif eye_Loc==5
        loc5(loc5_count) = iTrial;
        loc5_count = loc5_count + 1;
    elseif eye_Loc==6
        loc6(loc6_count) = iTrial;
        loc6_count = loc6_count + 1;
    elseif eye_Loc==7
        loc7(loc7_count) = iTrial;
        loc7_count = loc7_count + 1;
    elseif eye_Loc==8
        loc8(loc8_count) = iTrial;
        loc8_count = loc8_count + 1;
    end
end

%index loc_eye data
for i=1:8
    if i==1
        locXYstruct(i).locX = eye_catX(loc1,1:10:end);
        locXYstruct(i).locY = eye_catY(loc1,1:10:end);
    elseif i==2
        locXYstruct(i).locX = eye_catX(loc2,1:10:end);
        locXYstruct(i).locY = eye_catY(loc2,1:10:end);
    elseif i==3
        locXYstruct(i).locX = eye_catX(loc3,1:10:end);
        locXYstruct(i).locY = eye_catY(loc3,1:10:end);
    elseif i==4
        locXYstruct(i).locX = eye_catX(loc4,1:10:end);
        locXYstruct(i).locY = eye_catY(loc4,1:10:end);
    elseif i==5
        locXYstruct(i).locX = eye_catX(loc5,1:10:end);
        locXYstruct(i).locY = eye_catY(loc5,1:10:end);
    elseif i==6
        locXYstruct(i).locX = eye_catX(loc6,1:10:end);
        locXYstruct(i).locY = eye_catY(loc6,1:10:end);
    elseif i==7
        locXYstruct(i).locX = eye_catX(loc7,1:10:end);
        locXYstruct(i).locY = eye_catY(loc7,1:10:end);
    elseif i==8
        locXYstruct(i).locX = eye_catX(loc8,1:10:end);
        locXYstruct(i).locY = eye_catY(loc8,1:10:end);
    end
end

close all

plotDir = 'D:\\WTF_EYE\\Plots';

locXcat1 = [locXYstruct(1).locX;locXYstruct(5).locX ];
locXcat2 = [locXYstruct(2).locX;locXYstruct(6).locX ];
locXcat3 = [locXYstruct(3).locX;locXYstruct(7).locX ];
locXcat4 = [locXYstruct(4).locX;locXYstruct(8).locX ];

locYcat1 = [locXYstruct(1).locY;locXYstruct(5).locY ];
locYcat2 = [locXYstruct(2).locY;locXYstruct(6).locY ];
locYcat3 = [locXYstruct(3).locY;locXYstruct(7).locY ];
locYcat4 = [locXYstruct(4).locY;locXYstruct(8).locY ];


h=figure('units','normalized','outerposition',[0 0 1 1]);

for p=1:4
    
    subplot(1,4,p);
    
    if      p==1; theseDataX = locXcat1; theseDataY = locYcat1; 
    elseif  p==2; theseDataX = locXcat2; theseDataY = locYcat2;
    elseif  p==3; theseDataX = locXcat3; theseDataY = locYcat3;
    elseif  p==4; theseDataX = locXcat4; theseDataY = locYcat4;
    end
    
    for j=1:size(theseDataX,1)
        for i=1:size(theseDataY,2)
            plot(theseDataX(j,i), theseDataY(j,i),'linestyle','none','marker','.','MarkerSize',10); hold on
            set(gca,'xlim',[0,1024],'ylim',[0,768])
        end
    end
    
    pbaspect([1,1,1])
    
     
end

saveas(h,[plotDir '/' 'Eye_Movements_Loc_sj11.jpeg'],'jpeg')
