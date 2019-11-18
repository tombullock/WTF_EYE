%{
Author: Tom Bullock
Date: 11.14.19

List bad channels for specific subjects
%}

function badChannels = badChannelInfo(sjNum)

if sjNum==5
    badChannels = {'FC5','Iz','FC4','P9','TP8'};
elseif sjNum==999
    badChannels = {'P1','F3'};
else
    badChannels = {};
end

return
