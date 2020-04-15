%{
Beh_Rename_Files

Author: Tom Bullock, UCSB Attention Lab
Date: 10.30.19

Remove the randomized 4 digit code from the end of each raw beh file to make the files easier to work with.
Skip files that have already been renamed.
%}

function Beh_Rename_Files(behDir)

cd(behDir)

d=dir('sj*');

for i=1:length(d)
    if length(d(i).name)>28
        movefile([behDir '/' d(i).name],[behDir '/' d(i).name(1:24) '.mat'])
    end
end

cd ..

return