%{
IEM_Settings
Author: Tom Bullock (borrowed heavily from Joshua Foster)
Date: 11.16.19

%}

function em = IEM_Settings

% parameters to set
em.nChans = 8; % # of channels
em.nBins = em.nChans; % # of stimulus bins
em.nIter = 10; % # of iterations [UP THIS TO 10 EVENTUALLY]
em.nPerms = 3; % # of permutations (for permuted datasets) [KEEP THIS AT 3 OR FILES GET TOO BIG]
em.nBlocks = 3; % # of blocks for cross-validation
em.frequencies = [8,12;4,7]; % frequency bands to analyze
em.bands = {'Alpha','Theta'};
em.Fs = 256;
em.window = 4;
em.time = -.5*1000:1000/em.Fs:1.9961*1000; %   -500:4:2000; % time points of interest

% Specify basis set
nBins = em.nBins;
nChans = em.nChans;
em.sinPower = 7;
em.x = linspace(0, 2*pi-2*pi/nBins, nBins);
em.cCenters = linspace(0, 2*pi-2*pi/nChans, nChans);
em.cCenters = rad2deg(em.cCenters);
pred = sin(0.5*em.x).^em.sinPower; % hypothetical channel responses
pred = wshift('1D',pred,5); % shift the initial basis function
basisSet = nan(nChans,nBins);
for c = 1:nChans;
    basisSet(c,:) = wshift('1D',pred,-c); % generate circularly shifted basis functions
end
em.basisSet = basisSet; % save basis set to data structure

return