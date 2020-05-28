%==========================================================================
%{
calculateSlopes_Single_Freqs
Purpose: calculate slope values for single freq bandpass IEM analyses

Original Author:
Joshua J. Foster
joshua.james.foster@gmail.com
University of Chicago
August 17, 2015

Modified by Tom Bullock
UCSB Attention Lab
%}
%==========================================================================

clear all
close all

% setup directories
%root = pwd; out = 'AnalysisScripts/trunk/MATLAB';
%dRoot = [root(1:end-length(out)),'Data/'];
dRoot = '/home/bullock/WTF/Data/';

% which subjects?
subjects = [101:104 201:204 301:304 401:404 501:504 601:604 701:704 801:804 ...
    901:904 1001:1004 1301:1304 1601:1604 1901:1904 2001:2004 2101:2104 2201:2204 2301:2304 2401:2404 2501:2504];

% subjects = [3102,3105,3202,3205,3302,3305, 3402,3405, 3502,3505,3602,3605,3702,3705,...
%     3802,3805,3902,3905, 4002,4005,4102,4105,4202,4205]; 


for bandpassLoop=1
    
    if bandpassLoop==1
        name = '_SpatialTF_ALPHA.mat'; % name of files to be saved
        thisFreq = {'Alpha'};
        freqBandpass = [8 12];
        saveName = 'CTFslopes_Single_Freq_ALPHA.mat';
    else
        name = '_SpatialTF_THETA.mat';
        thisFreq = {'Theta'};
        freqBandpass = [4 7];
        saveName = 'CTFslopes_Single_Freq_THETA.mat';
    end
    
    matlabpool open 72
    parfor iSub=1:length(subjects)
        
        % clear stuff before next loop
        em = []; pDat = []; pSl = []; rDat = []; rSl = []; d= []; dat = [];rPw=[]; pPw=[];
        
        sn = subjects(iSub);
        
        % grab subject's data (*_SpatialTF_Permed also contains real DTFs)
        fName = [dRoot,num2str(sn), name];
        tmp = load(fName);
        em = tmp.em;
        tmp.em = [];
        rDat.evoked = em.tfs.evoked;
        rDat.total = em.tfs.total;  
        %     rDat.evoked = squeeze(mean(mean(em.tfs.evoked,2),4)); % average across iteration and cross-validation blocks rDAT = REAL DATA
        %     rDat.total = squeeze(mean(mean(em.tfs.total,2),4));
        pDat.evoked = squeeze(mean(em.permtfs.evoked,2)); % average across iterations pDAT = PERMUTED DATA
        pDat.total = squeeze(mean(em.permtfs.total,2));
        
        % Specify properties
        nChans = em.nChans; % # of location channels
        nPerms = 10; % % of permutations
        %nSamps = em.dSamps; % # of samples (after downsampling)
        %nFreqs = em.nFreqs; % # of frequencies
        nSamps = size(rDat.total,1);
        nFreqs = 1;
        
        % SPECIFY X-values (TOM ADDITION)
        %thisX = 1:5; % foster uses this
        thisX = 0:45:180; % WTF use real angular values
        
        % real evoked data
        for f = 1:nFreqs
            for samp = 1:nSamps;
                dat = rDat.evoked(samp,:);
                %dat = squeeze(rDat.evoked(f,samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                rSl.evoked(f,samp)= fit(1);
                rPw.evoked(f,samp)=dat(5);
            end
        end
        
        % real total data
        for f = 1:nFreqs
            for samp = 1:nSamps;
                dat = rDat.total(samp,:);
                %dat = squeeze(rDat.evoked(f,samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                rSl.total(f,samp)= fit(1);
                rPw.total(f,samp)=dat(5);
            end
        end
        
        
        % permuted evoked data
        for perm = 1:nPerms
            disp(['Evoked Perm: ' num2str(perm)])
            for f = 1:nFreqs
                for samp = 1:nSamps;
                    dat = squeeze(pDat.evoked(perm,samp,:));
                    %dat = squeeze(pDat.evoked(f,perm,samp,:));
                    x = thisX;
                    d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                    fit = polyfit(x,d,1);
                    pSl.evoked(f,samp,perm)= fit(1) ;
                    pPw.evoked(f,samp,perm)=dat(5);
                end
            end
        end
        
        % permuted total data
        for perm = 1:nPerms
            disp(['Total Perm: ' num2str(perm)])
            for f = 1:nFreqs
                for samp = 1:nSamps;
                    dat = squeeze(pDat.total(perm,samp,:));
                    %dat = squeeze(pDat.total(f,perm,samp,:));
                    x = thisX;
                    d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                    fit = polyfit(x,d,1);
                    pSl.total(f,samp,perm)= fit(1) ;
                    pPw.total(f,samp,perm)=dat(5);
                end
            end
        end
        
        % save slope matrices
        filename = [dRoot,num2str(sn),saveName];
        parsave(filename,rSl,pSl,rPw,pPw);
        
    end
    
    matlabpool close
    
end

clear all
close all


