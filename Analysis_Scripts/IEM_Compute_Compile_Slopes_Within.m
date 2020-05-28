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

UPDATED AGAIN 05.02.20 BY TOM
%}

%function IEM_Compute_Compile_Slopes_Within(subjects)

%if nargin==0
subjects = [1:7,9:14,16:20,22:27,31];
%end

%clear
%close all

% set dirs
rDir = '/home/waldrop/Desktop/WTF_EYE';
sourceDir = [rDir '/' 'IEM_Results_TT_Within_Tmp' ];
%destDir = [rDir '\\' 'IEM_Slopes_TT_Within'];
destDirCompiled = [rDir '/' 'Data_Compiled'];

% subjects 
%subjects = [4,5];

% specify x values
thisX = 0:45:180; 

% sub loop
for iSub=1:length(subjects)
    sjNum = subjects(iSub)
    
    for iCond=1:4
       
        % load the "within" data only
        load([sourceDir '/' sprintf('sj%02d_cond%02d_IEM.mat',sjNum,iCond)],'em_within')
        
        % real total data
        for f = 1:size(em_within.tfs.total,1)
            for samp = 1:size(em_within.tfs.total,2)
                dat = squeeze(em_within.tfs.total(f,samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                %em_within.rSl.total(f,samp)= fit(1);
                %allSlopes_real_total(iSub,iCond,f,samp) = fit(1);
                rSl_total(iSub,iCond,f,samp) = fit(1); 
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
                %allSlopes_real_evoked(iSub,iCond,f,samp) = fit(1);
                rSl_evoked(iSub,iCond,f,samp) = fit(1);
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
                %allSlopes_perm_total(iSub,iCond,f,samp) = fit(1);
                pSl_total(iSub,iCond,f,samp) = fit(1);
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
                %allSlopes_perm_evoked(iSub,iCond,f,samp) = fit(1);
                pSl_evoked(iSub,iCond,f,samp) = fit(1);
            end
        end
        
        clear em_within
        
    end
    
    
end




% convert into R "long" format

% freqs loop (alpha, theta)
for d=1:4
    for f=1:2
        for t=1:640
            dumMem = [zeros(25,1), ones(25,1), zeros(25,1), ones(25,1)];
            dumEyes = [zeros(25,1), zeros(25,1), ones(25,1), ones(25,1)];
            
            % create column vector of subjects
            sjNums = 1:size(rSl_total,1);
            sjNums = [sjNums, sjNums, sjNums, sjNums]';
            
            % isolate data
            if      d==1; theseData = rSl_total(:,:,f,t);
            elseif  d==2; theseData = rSl_evoked(:,:,f,t);
            elseif  d==3; theseData = pSl_total(:,:,f,t);
            elseif  d==4; theseData = pSl_evoked(:,:,f,t);
            end

            % write out data
            if      d==1; rSl_total_LF(:,:,f,t) = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
            elseif  d==2; rSl_evoked_LF(:,:,f,t) = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
            elseif  d==3; pSl_total_LF(:,:,f,t) = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
            elseif  d==4; pSl_evoked_LF(:,:,f,t) = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
            end
            
        end
    end
end




% save all long and regular format data
save([destDirCompiled '/' 'IEM_Slopes_Within.mat'],...
    'rSl_evoked',...
    'rSl_total',...
    'pSl_evoked',...
    'pSl_total',...
    'rSl_evoked_LF',...
    'rSl_total_LF',...
    'pSl_evoked_LF',...
    'pSl_total_LF','-v7.3')
        
  
%return
