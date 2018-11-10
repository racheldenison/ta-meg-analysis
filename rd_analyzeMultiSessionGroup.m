% rd_analyzeMultiSessionGroup.m

%% setup
exptDir = '/Local/Users/denison/Data/TANoise/MEG';

t = -1500:5700;
tidx = 1:6701;

sessionDirsAll = {{'R0817_20171212','R0817_20171213'},...
                {'R1187_20180105','R1187_20180108'},...
                {'R0983_20180111','R0983_20180112'},...
                {'R0898_20180112','R0898_20180116'},...
                {'R1021_20180208','R1021_20180212'},...
                {'R1103_20180213','R1103_20180215'},...
                {'R0959_20180219','R0959_20180306'}};

subjects = 1:numel(sessionDirsAll);
nSubjects = numel(subjects);

analysisDir = sprintf('%s/Group/mat', exptDir);
peaksFileName = sprintf('%s/gN%d_peaks.mat', analysisDir, nSubjects);
itpcFileName = sprintf('%s/gN%d_itpcAtt.mat', analysisDir, nSubjects);

saveAnalysis = 0;

%% get data
data = [];
for iSubject = 1:nSubjects
    subject = subjects(iSubject);
    sessionDirs = sessionDirsAll{subject};
    
    [data(iSubject).itpc, data(iSubject).peaks] = rd_plotMultiSessionData(sessionDirs);
    
    close all
end

%% organize data
itpc = [];
for iSubject = 1:nSubjects
    itpc(:,:,iSubject) = data(iSubject).itpc(tidx,:);
end

%% save analysis
if saveAnalysis
    save(peaksFileName, 'exptDir', 'sessionDirsAll', 'subjects', 'peaks')
    save(itpcFileName, 'exptDir', 'sessionDirsAll', 'subjects', 't', 'itpc')
end