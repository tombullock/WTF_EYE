%{
Correlate_Beh_With_CTFs
Author: Tom Bullock
Date: 05.19.20

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';
destDir = sourceDir;

% load slope data
load([sourceDir '/' 'IEM_Slopes_Within.mat'],'rSl_total')

% load beh data
load([sourceDir '/' 'Modelling_Data.mat'])

% remove bad subs
badSubIdx = 11;
modelSD(badSubIdx,:) = [];
modelGuess(badSubIdx,:) = [];


% create times vector for reference
times = linspace(-500,2000,640);

% isolate slopes times to average over (change these values if needed)
peakTimes = [200,400]; % in milliseconds!

[~,peakTimeIdx(1)] = min(abs(peakTimes(1)-times));
[~,peakTimeIdx(2)] = min(abs(peakTimes(2)-times));

% get data
thesePeaks = squeeze(mean(rSl_total(:,:,1,peakTimeIdx(1):peakTimeIdx(2)),4));

% remove bad subs
thesePeaks(badSubIdx,:) = [];

%plot(mean(thesePeaks))



% plot correlations between slope and behavioral measures
h=figure;

cnt=0;
for d = 1:2
    
    if d == 1; theseBehData = modelSD;
    else; theseBehData = modelGuess;
    end
    
    % if calculating delta
    
    for p=1:4
        
        cnt=cnt+1;
        
        subplot(2,4,cnt);
        scatter(thesePeaks(:,p),theseBehData(:,p));
        lsline
        
        pbaspect([1,1,1])
        
        [rho,pval] = corr(thesePeaks(:,p),theseBehData(:,p));
        
        title(['Rho = ' num2str(rho) '   pVal = ' num2str(pval)])
        
    end
end

% plot correlations between slope delta (spatial fixed - spatial moved) and
% behavioral measures

% if calculating slope delta, keep this line
thesePeaksDeltaSpatial = thesePeaks(:,1) - thesePeaks(:,3);
thesePeaksDeltaColor = thesePeaks(:,2) - thesePeaks(:,4);

sdSpatial = modelSD(:,1) - modelSD(:,3);
sdColor = modelSD(:,2) - modelSD(:,4);


h=figure;
subplot(1,2,1)
scatter(thesePeaksDeltaSpatial,sdSpatial)
lsline
pbaspect
[rho,pval] = corr(thesePeaksDeltaSpatial,sdSpatial)     
title(['Rho = ' num2str(rho) '   pVal = ' num2str(pval)])

subplot(1,2,2)
scatter(thesePeaksDeltaColor,sdColor)
lsline
pbaspect
[rho,pval] = corr(thesePeaksDeltaColor,sdColor)     
title(['Rho = ' num2str(rho) '   pVal = ' num2str(pval)])



