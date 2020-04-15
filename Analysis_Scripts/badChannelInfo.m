%{
Author: Tom Bullock
Date: 11.14.19

List bad channels for specific subjects
%}

function badChannels = badChannelInfo(sjNum)

if sjNum==1
    badChannels = {'O1','T7','T8','C6'};
elseif sjNum==5
    badChannels = {'Iz','P9','FC5','PO7','FC4','P10'};    
elseif sjNum==6
    badChannels = {'F7','Fpz','P5'};
elseif sjNum==2
    badChannels = {'AF7','Fp1'};
elseif sjNum==9
    badChannels = {'C5','C3','CP3','Oz'};
elseif sjNum==10
    badChannels = {'F4','TP7','Fp2','PO8','TP8','P10'};
elseif sjNum==12
    badChannels = {'FC5','P10'};
elseif sjNum==13
    badChannels = {'PO4','O1'};
elseif sjNum==14
    badChannels = {'AF8','F4'};
elseif sjNum==16
    badChannels = {'FT8'};
elseif sjNum==17
    badChannels = {'CPz'};
elseif sjNum==18
    badChannels = {'CP6','P9'};
elseif sjNum==19
    badChannels = {'F8'};
elseif sjNum==22
    badChannels = {'TP7'};
elseif sjNum==25
    badChannels = {'P2','F7'};
elseif sjNum==31
    badChannels = {'Fp2'};
else
    badChannels = {};
end

return
