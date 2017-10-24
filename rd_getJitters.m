function behav = rd_getJitters(exptDir, sessionDir, runs, behav)

%% Setup
if nargin==0
    exptDir = '/Local/Users/denison/Data/TAContrast/MEG';
    sessionDir = 'R1187_20170510';
    runs = 21:23;
    behav = [];
    excludeBlanks = 0;
end

behavDir = sprintf('%s/Behavior/%s/data', exptDir(1:end-4), sessionDir);

nRuns = numel(runs);

%% Get data
for iRun = 1:nRuns
    run = runs(iRun);
    behavFile = dir(sprintf('%s/*%d.mat', behavDir, run));
    b(iRun) = load(sprintf('%s/%s', behavDir, behavFile.name));
end

%% Get jitters
jitSeq = [];
for iRun = 1:nRuns
    jitSeq = [jitSeq b(iRun).stimulus.jitSeq];
end

%% Exclude blanks
% if excludeBlanks
%     cueCondIdx = strcmp(behav.responseData_labels, 'cue condition');
%     cueCond = behav.responseData_all(:,cueCondIdx);
%     wBlank = cueCond==1;
%     jitSeq(wBlank) = [];
% end
% 
% %% Output
% behav.jitSeq = jitSeq';

%% Add to responseData
behav.responseData_all(:,end+1) = jitSeq';
behav.responseData_labels{end+1} = 'jitter';