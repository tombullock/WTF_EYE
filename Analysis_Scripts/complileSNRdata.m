subjects = [1:8,10:35];
conds = 1:2;

root = '/home/garrett/WTF_Bike/';
SNR_Dir = [root 'SNR/'];
outDir = [root 'Data_Compiled/'];

SNR = [];
Exercise = [];
electrode_snr = [];
SubID = [];
Time = [];
for iSub = 1:length(subjects)
    sjNum = subjects(iSub);
    SubID = [SubID;repmat(sjNum,6,1)];
    for c = conds
        noA = load([SNR_Dir sprintf('sj%02d_exerCon%02d_SNR_Hilbert_noAlpha.mat',sjNum,c)],'rawHilbert'); %% This is a electrode X timepoints matrix of non-alpha power
        A = load([SNR_Dir sprintf('sj%02d_exerCon%02d_SNR_Hilbert_Alpha.mat',sjNum,c)],'avgPower'); %% This is a electrode X timepoints matrix of alpha power
        A = A.avgPower; 
        
        go = ' ';
        fo = squeeze(mean(A,1));
        foo(c,iSub,:) = fo - mean(fo(1:125));
        
        if go == ' '
            continue
        end
        
        noA = noA.rawHilbert; %really the std of the rawHilbert
        
        time = linspace(-500,2000,size(A,2));
        for t = 1:3
            if t == 1
                tois = [find(abs(time-0)==min(abs(time-0))):find(abs(time-250)==min(abs(time-250)))]; % encoding
            elseif t== 2
                tois = [find(abs(time-250)==min(abs(time-250))):find(abs(time-2000)==min(abs(time-2000)))]; % delay period
            elseif t == 3
                tois = [find(abs(time-0)==min(abs(time-0))):find(abs(time-2000)==min(abs(time-2000)))]; % whole period
            end
            
            electrode_snr(iSub,:,:) = A./noA; % save for looking at spatial SNR
            temp_snr = squeeze(mean(electrode_snr(iSub,:,:))); %% This is basically it. 
           
            temp_snr = squeeze(mean(temp_snr(tois))); %% I did two different time windows for this, not necessary
            SNR = [SNR;temp_snr];
            Time = [Time;t];
            if c== 1
                Exercise = [Exercise;0];
            elseif c == 2
                Exercise = [Exercise;1];
            end
        end
        
        save([outDir sprintf('exerCon%02d_SNR_all.mat',c)],'electrode_snr');
    end
end

data = table(SubID,SNR,Exercise,Time);
writetable(data,[SNR_Dir 'SNR_data_table.csv'],'Delimiter',',')