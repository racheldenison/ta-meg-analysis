function rd_writeTAMEGDataFile(dataFileName)

% dataFileName = '/Local/Users/denison/Data/TADetectDiscrim/MEG/Group/mat/ebi_ft/gN16_singleTrials_30Hz_topChannels5_allTrials_normStim_itpc-single.mat';

% print all the file names
% dataDir = '/Local/Users/denison/Data/TADetectDiscrim/MEG/Group/mat/ebi_ft'
% dir(sprintf('%s/gN16*PAAUT*.txt', dataDir))

%% Load data
D = load(dataFileName);
info = D.info;
data = D.data;

subjects = info.subjects;
measure = info.measure;
eventTimes = info.eventTimes;
try
    condNames = info.condNames;
catch
    condNames = info.trigNames;
end

nSubjects = numel(subjects);

if eventTimes==0
    paaut = true;
else
    paaut = false;
end

%% Setup
if paaut
    twin = [info.t(1) info.t(end)];
else
    twin = eventTimes([2 5]); % [500 3100];
end

condIdx = 1:8; % do not include blank, so we have 2x2x2
condStr = '';
for iCond = condIdx
    condStr = sprintf('%s %s', condStr, condNames{iCond});
end

%% Determine time points
switch measure
    case {'tf','stf','stf-single'}
        if paaut
            t = info.t;
        else
            t = round(info.stfToi*1000);
        end
    otherwise
        t = info.t;
end
tidx = find(t==twin(1)):find(t==twin(2));
nT = numel(tidx);

%% Select data if needed
extraStr = '';
switch measure
    case {'ts','ts-single'}
        if ~paaut
            data = squeeze(data(:,1,:,:)); % top channel only
        end
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
fprintf(fileID, '%s %s %s\n','subject','time', condStr(2:end));

for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    
    for iT = 1:nT
        fprintf(fileID,'%s %d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f\n', ...
            subject, t(tidx(iT)), data(tidx(iT),condIdx,iSubject));
    end
end
    
fclose(fileID);