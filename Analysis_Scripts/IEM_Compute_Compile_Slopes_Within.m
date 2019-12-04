%{
IEM_Compute_Slopes_Within
Purpose: calculate slope values for single freq bandpass IEM analyses

Original Author:
Joshua J. Foster
joshua.james.foster@gmail.com
University of Chicago
August 17, 2015

Modified by Tom Bullock
UCSB Attention Lab
%}

function IEM_Compute_Compile_Slopes_Within(subjects)

%clear
%close all

% set dirs
rDir = 'C:\\Users\\BOSS-EEG\\Desktop\\WTF_EYE';
sourceDir = [rDir '\\' 'IEM_Results_TT_Within' ];
%destDir = [rDir '\\' 'IEM_Slopes_TT_Within'];
destDirCompiled = [rDir '\\' 'Data_Compiled'];

% subjects 
%subjects = [4,5];

% specify x values
thisX = 0:45:180; 

% sub loop
for iSub=1:length(subjects)
    sjNum = subjects(iSub);
    
    for iCond=1:4
       
        % load the "within" data only
        load([sourceDir '\\' sprintf('sj%02d_cond%02d_IEM.mat',sjNum,iCond)],'em_within')
        
        % real total data
        for f = 1:size(em_within.tfs.total,1)
            for samp = 1:size(em_within.tfs.total,2)
                dat = squeeze(em_within.tfs.total(f,samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                %em_within.rSl.total(f,samp)= fit(1);
                allSlopes_real_total(iSub,iCond,f,samp) = fit(1);
            end
        end
        
        % real evoked data
        for f = 1:size(em_within.tfs.evoked,1)
            for samp = 1:size(em_within.tfs.evoked,2)
                dat = squeeze(em_within.tfs.evoked(f,samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                %em_within.rSl.evoked(f,samp)= fit(1);
                allSlopes_real_evoked(iSub,iCond,f,samp) = fit(1);
            end
        end
        

        % perm total data
        for f = 1:size(em_within.tfs_perm.total,1)
            for samp = 1:size(em_within.tfs_perm.total,2)
                dat = squeeze(em_within.tfs_perm.total(f,samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                %em_within.pSl.total(f,samp)= fit(1);
                allSlopes_perm_total(iSub,iCond,f,samp) = fit(1);

            end
        end
        
        % perm evoked data
        for f = 1:size(em_within.tfs_perm.evoked,1)
            for samp = 1:size(em_within.tfs_perm.evoked,2)
                dat = squeeze(em_within.tfs_perm.evoked(f,samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                %em_within.pSl.evoked(f,samp)= fit(1);
                allSlopes_perm_evoked(iSub,iCond,f,samp) = fit(1);
            end
        end
        
        clear em_within
        
    end
    
end


save([destDirCompiled '\\' 'IEM_Slopes_Within.mat'],...
    'allSlopes_perm_evoked',...
    'allSlopes_perm_total',...
    'allSlopes_real_evoked',...
    'allSlopes_real_total')

return
