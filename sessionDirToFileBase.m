function fileBase = sessionDirToFileBase(sessionDir)

% function fileBase = sessionDirToFileBase(sessionDir)
%
% example:
% sessionDir = 'R0817_20150504';
% fileBase = 'R0817_TADeDi_5.4.15';

s = strsplit(sessionDir,'_');

rnum = s{1};
date = s{2};

year = date(3:4);
month = date(5:6);
day = date(7:8);

if month(1)=='0'
    month = month(2);
end
if day(1)=='0'
    day = day(2);
end

fileBase = sprintf('%s_TADeDi_%s.%s.%s', rnum, month, day, year);