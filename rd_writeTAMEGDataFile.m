function rd_writeTAMEGDataFile(dataFileName)

% dataFileName = '/Local/Users/denison/Data/TADetectDiscrim/MEG/Group/mat/ebi_ft/gN16_singleTrials_30Hz_topChannels5_allTrials_normStim_itpc-single.mat';

% print all the file names
% dataDir = '/Local/Users/denison/Data/TADetectDiscrim/MEG/Group/mat/ebi_ft'
% dir(sprintf('%s/gN16*.txt', dataDir))

%% Load data
D = load(dataFileName);
info = D.info;
data = D.data;

subjects = info.subjects;
measure = info.measure;
eventTimes = info.eventTimes;
trigNames = info.trigNames;

nSubjects = numel(subjects);

%% Setup
condIdx = 1:8; % do not include blank, so we have 2x2x2
twin = eventTimes([2 5]); % [500 3100];

trigStr = '';
for iCond = condIdx
    trigStr = sprintf('%s %s', trigStr, trigNames{iCond});
end

%% Determine time points
switch measure
    case {'tf','stf','stf-single'}
        t = round(info.stfToi*1000);
    otherwise
        t = info.t;
end
tidx = find(t==twin(1)):find(t==twin(2));
nT = numel(tidx);

%% Select data if needed
extraStr = '';
switch measure
    case {'ts','ts-single'}
        data = squeeze(data(:,1,:,:)); % top channel only
    case {'tf','stf','stf-single'}
        f = info.stfFoi;
        fwin = [8 12]; % select alpha range
        fidx = find(f==fwin(1)):find(f==fwin(2));
        data = squeeze(nanmean(data(fidx,:,:,:)));
        extraStr = '8-12Hz_';
end

%% Write text file (wide format)
textFileName = sprintf('%s_%s%d-%dms.txt', dataFileName(1:end-4), extraStr, twin);

fileID = fopen(textFileName,'w');
fprintf(fileID, '%s %s %s\n','subject','time', trigStr(2:end));

for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    
    for iT = 1:nT
        fprintf(fileID,'%s %d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n', ...
            subject, t(tidx(iT)), data(tidx(iT),condIdx,iSubject));
    end
end
    
fclose(fileID);