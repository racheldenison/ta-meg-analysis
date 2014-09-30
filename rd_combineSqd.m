function outFile = rd_combineSqd(dataDir, outName, tag, delimiter)
%
% function outFile = rd_combineSqd([dataDir], [outName], [tag], [delimiter])
%
% Append multiple sqd files. Finds all the sqd files in a directory that
% have a given tag.
%
% INPUTS:
% dataDir is the full path to the directory with the sqds to combine.
%   Default is the current directory.
% outName is the name of the combined sqd file. This will be saved in
%   dataDir. Default is 'combined_[current time].sqd'.
% tag is an extension that follows a delimiter and appears just before .sqd
%   in the names of the files you want to combine. E.g. 'run' or 'edts'.
%   Default is '', which means combine all the files in the directory.
% delimiter is a string that separates the tag from the previous part of
%   the file name. Default is '_'.
%
% Rachel Denison
% September 2014

%% deal with inputs
if nargin < 4
    delimiter = '_';
end
if nargin < 3
    tag = ''; % take all the files
end
if nargin < 2
    outName = sprintf('combined_%s.sqd', datestr(now,30));
end
if nargin < 1
    dataDir = pwd;
end

%% setup
% this "tag" is a little bit of a complicated idea
% it means that the string after the last delimiter should include this
% string
% delimiter = '_';
% tag = 'run';
% 
% dataDir = '/Local/Users/denison/Data/TAPilot/MEG/R0890_20140806/Runs';
% outName = 'R0890_TAPilot_8.06.14_run1and2.sqd';
outFile = sprintf('%s/%s', dataDir, outName);

%% get all sqd file names from the directory
sqdFiles = dir([dataDir '/*.sqd']);

%% find the files with the right tags
sqdToCombine = {};
for iFile = 1:numel(sqdFiles)
    sqdName = sqdFiles(iFile).name;
    tagStr = rd_getTag(sqdName, delimiter);
    if ~isempty(strfind(tagStr,tag)) || strcmp(tag,'')
        sqdToCombine = [sqdToCombine, sqdName]; 
    end
end

disp(sqdToCombine')

%% check with the user that this looks good
go = input('List of sqd files to combine look ok? [y,n]','s');
if ~strcmp(go,'y')
    fprintf('\nNot combining files ... exiting\n\n')
    return
end

%% make combined sqd file
sourceFile = sprintf('%s/%s', dataDir, sqdToCombine{1});
for iFile = 1:numel(sqdToCombine)
    sqdNow = sprintf('%s/%s', dataDir, sqdToCombine{iFile});
    data = sqdread(sqdNow);
    sqdwrite(sourceFile, outFile, 'action', 'append', 'data', data);
end
