%{
IEM_Settings
Author: Tom Bullock (borrowed heavily from Joshua Foster)
Date created: 11.16.19
Date updated: 10.23.20
Inputs: none

Outputs: em (encoding model) structure containing all IEM settings

%}

function em = IEM_Fixed_Model_LOA_Settings

% parameters to set
em.nChans = 8; % # of hypothetical orientation/location channels
em.nBins = em.nChans; % # of stimulus bins (typically equal to nChans)
em.nIter = 1; % # of iterations
%em.nPerms = 5; % # of permutations 
em.nBlocks = 3; % # of blocks for cross-validation
em.frequencies = [8,12]; % frequency bands to analyze
em.bands = {'Alpha'}; % frequency band labels
em.Fs = 256; % sample rate
em.window = 4; % ?remove
em.time = -.5*1000:1000/em.Fs:1.9961*1000; % time window e.g. [-500:4:2000]

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