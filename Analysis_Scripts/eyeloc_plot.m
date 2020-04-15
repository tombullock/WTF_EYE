close all

for iTrial=1:length(trial_data)
    eye_Loc = trial_data(iTrial).eyeLoc;
    if eye_Loc==1
        plot (trial_data(iTrial).eyeLocX,trial_data(iTrial).eyeLocY,'linestyle','none','marker','.','MarkerSize',10,'color','r'); hold on
    elseif eye_Loc==2
        plot (trial_data(iTrial).eyeLocX,trial_data(iTrial).eyeLocY,'linestyle','none','marker','.','MarkerSize',10,'color',[255,140,0]./255)
    elseif eye_Loc==3
        plot (trial_data(iTrial).eyeLocX,trial_data(iTrial).eyeLocY,'linestyle','none','marker','.','MarkerSize',10,'color','y');
    elseif eye_Loc==4
        plot (trial_data(iTrial).eyeLocX,trial_data(iTrial).eyeLocY,'linestyle','none','marker','.','MarkerSize',10,'color','g');
    elseif eye_Loc==5
        plot (trial_data(iTrial).eyeLocX,trial_data(iTrial).eyeLocY,'linestyle','none','marker','.','MarkerSize',10,'color','b');
    elseif eye_Loc==6
        plot (trial_data(iTrial).eyeLocX,trial_data(iTrial).eyeLocY,'linestyle','none','marker','.','MarkerSize',10,'color',[148,0,211]./255);
    elseif eye_Loc==7
        plot (trial_data(iTrial).eyeLocX,trial_data(iTrial).eyeLocY,'linestyle','none','marker','.','MarkerSize',10,'color',[255,20,147]./255);
    elseif eye_Loc==8
        plot (trial_data(iTrial).eyeLocX,trial_data(iTrial).eyeLocY,'linestyle','none','marker','.','MarkerSize',10,'color',[128,128,128]./255);
    end
    set(gca,'xlim',[-300,300],'ylim',[-300,300])
end