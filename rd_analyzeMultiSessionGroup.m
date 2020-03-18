% rd_analyzeMultiSessionGroup.m

%% setup
% exptDir = '/Local/Users/denison/Data/TANoise/MEG';
exptDir = pathToTANoise('MEG');

t = -1500:5700;
tidx = 1:6701;

collapseSessions = 0;

sessionDirsAll = {{'R0817_20171212','R0817_20171213'},...
                {'R1187_20180105','R1187_20180108'},...
                {'R0983_20180111','R0983_20180112'},...
                {'R0898_20180112','R0898_20180116'},...
                {'R1021_20180208','R1021_20180212'},...
                {'R1103_20180213','R1103_20180215'},...
                {'R0959_20180219','R0959_20180306'},...
                {'R1373_20190723','R1373_20190725'},...
                {'R1452_20190717','R1452_20190718'},...
                {'R1507_20190702','R1507_20190705'}};

subjects = 1:numel(sessionDirsAll);
nSubjects = numel(subjects);

if collapseSessions
    sessStr = '';
else
    sessStr = '_bySession';
end

analysisDir = sprintf('%s/Group/mat', exptDir);
peaksFileName = sprintf('%s/gN%d_peaks%s.mat', analysisDir, nSubjects, sessStr); % from wSpecAll
itpcFileName = sprintf('%s/gN%d_itpcAtt%s.mat', analysisDir, nSubjects, sessStr); % from wSpecAtt
tsrFileName = sprintf('%s/gN%d_tsrAtt%s.mat', analysisDir, nSubjects); % time segment reliability

saveAnalysis = 0;

%% itpc and peaks
% get data
data = [];
for iSubject = 1:nSubjects
    subject = subjects(iSubject);
    sessionDirs = sessionDirsAll{subject}
    
    [data(iSubject).itpc, data(iSubject).peaks] = rd_plotMultiSessionData(sessionDirs, collapseSessions);
    
    close all
end

% organize data
itpc = [];
if collapseSessions
    itpc_dims = {'time','att','subject'};
    for iSubject = 1:nSubjects
        itpc(:,:,iSubject) = data(iSubject).itpc(tidx,:);
    end
else
    itpc_dims = {'time','att','subject','session'};
    for iSubject = 1:nSubjects
        itpc(:,:,iSubject,:) = data(iSubject).itpc(tidx,:,:);
    end
end

for iSubject = 1:nSubjects
    peaks(iSubject) = data(iSubject).peaks;
end

%% time segment reliability
tsrType = 'dist'; % 'corr','dist'

% get data
data = [];
for iSubject = 1:nSubjects
    subject = subjects(iSubject);
    sessionDirs = sessionDirsAll{subject};
    
    for iSession = 1:numel(sessionDirs)
        sessionDir = sessionDirs{iSession};
        switch tsrType
            case 'corr'
                data(iSubject,iSession).tsr = rd_timeSegmentReliabilityCorr(sessionDir);
            case 'dist'
                data(iSubject,iSession).tsr = rd_timeSegmentReliabilityDist(sessionDir);
        end
    end
end

% organize data
nCue = numel(data(1,1).tsr);

tsr = [];
for iSubject = 1:nSubjects
    for iSession = 1:numel(sessionDirs)
        for iCue = 1:nCue
            switch tsrType
                case 'corr'
                    tsr(:,iCue,iSubject) = mean(mean(data(iSubject,iSession).tsr{iCue}(tidx,:,:),3),2);
                case 'dist'
                    tsr(:,iCue,iSubject) = mean(data(iSubject,iSession).tsr{iCue}(tidx,:,:),2);
            end
            
        end
    end
end

figure
hold on
shadedErrorBar(t(tidx),mean(tsr(:,1,:),3),std(tsr(:,1,:),0,3)/sqrt(nSubjects),'b',1)
shadedErrorBar(t(tidx),mean(tsr(:,2,:),3),std(tsr(:,2,:),0,3)/sqrt(nSubjects),'r',1)
xlabel('Time (ms)')
ylabel(sprintf('Time segment %s', tsrType))
eventTimes = [0 1050 1350 2300];
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k')
end

%% save analysis
if saveAnalysis
    save(peaksFileName, 'exptDir', 'sessionDirsAll', 'subjects', 'peaks')
    save(itpcFileName, 'exptDir', 'sessionDirsAll', 'subjects', 't', 'tidx', 'itpc','itpc_dims')
    save(tsrFileName, 'exptDir', 'sessionDirsAll', 'subjects', 't', 'tidx', 'tsr')
end
