%{
IEM_Cross_TT
Author: Tom Bullock (using code borrowed from Joshua Foster)
Date: 11.18.19

Train and test within and across all conditions

To Do: build in the permuted data analysis
%}

%function IEM_Cross_TT_Calculate_Slope(subjects)

%clear 
%close all

% set source dir 
sourceDir  = '/home/waldrop/Desktop/WTF_EYE/IEM_Results_TT_Cross';
destDir = sourceDir;

% select subjects
%subjects = [1:7,9:12,16:19,20,22,23]; %orig
subjects = [1:7,9:14,16:20,22:27,31];


% loop through all training\testing combinations
for iSub=1:length(subjects)
    
    sjNum=subjects(iSub)
    
    em=[]; tf_total=[]; tf_total_cross=[]; tf_total_cross_perm=[];
    
    load([sourceDir '/' sprintf('sj%02d_IEM_Alpha_CrossTT.mat',sjNum) ])
    
    tf_total_cross = em.tfs.total_real;
    tf_total_cross_perm = em.tfs.total_perm;
    
    em.rSl = [];
    em.pSl = [];

    % calculate and save slopes (real then perm)
    thisX = 0:45:180; % WTF use real angular values
    for permReal=1:2
        
        if permReal==1; theseData = tf_total_cross;
        elseif permReal==2; theseData = tf_total_cross_perm;
        end
        
        for iCross=1:16
            for trSamps=1:40
                for teSamps=1:40
                    %dat=rDat.total(iCross,trSamps,teSamps,:);
                    dat = theseData(iCross,trSamps,teSamps,:);
                    x = thisX;
                    d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                    
                    
                    %fit = trapz(d);
                    
                    fit = polyfit(x,d,1);
                    
                    if      permReal==1; em.rSl.total(iCross,trSamps,teSamps)=fit(1);
                    elseif  permReal==2; em.pSl.total(iCross,trSamps,teSamps)=fit(1);
                    end
                end
            end
        end
    end
    
%     % save em
%     em.tfs.total_real = tf_total_cross;
%     em.tfs.total_perm = tf_total_cross_perm;

    quant = 'Slope';
    
    save([destDir '/' sprintf('sj%02d_IEM_Alpha_CrossTT.mat',sjNum)],'em','quant');
    
end


%clear 
%close all

%return