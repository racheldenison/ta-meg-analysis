% rd_writeAttDataFile.m

%% i/o
exptDir = pathToTANoise('MEG');
analysisDir = sprintf('%s/Group/mat', exptDir);

% from rd_plotPeaks.m
analysisName = 'gN10_itpcAttPeakVals_20Hz_bySession';

textfileName = 'TANoise_N10_itpcAttPeakVals_20Hz';

%% load data
D = load(sprintf('%s/%s.mat', analysisDir, analysisName));

%% setup
data0 = D.peakDataAve; % [cueT1/T2 peak1/2 subject session]

% rearrange from cueT1/T2 to att/unatt
data(1,1,:,:) = data0(1,1,:,:); % T1 att
data(2,1,:,:) = data0(2,1,:,:); % T1 unatt
data(1,2,:,:) = data0(2,2,:,:); % T2 att
data(2,2,:,:) = data0(1,2,:,:); % T2 unatt

nAtt = size(data,1);
nTarget = size(data,2);
nSubject = size(data,3);
nSession = size(data,4);

targetNames = {'t1','t2'};
attNames = {'att','unatt'};

%% write text file
fileID = fopen(sprintf('%s/%s.txt', analysisDir, textfileName),'w');
fprintf(fileID,'%s %s %s %s %s\n','subject','session',...
    'target','att','measure');

for iSub = 1:nSubject
    subject = iSub;
    
    for iSess = 1:nSession
        session = iSess;
        
        for iT = 1:nTarget
            target = targetNames{iT};
            
            for iAtt = 1:nAtt
                att = attNames{iAtt};
                
                val = data(iAtt,iT,iSub,iSess);
                
                fprintf(fileID,'%d %d %s %s %1.4f\n', ...
                    subject, session, target, att, val);
            end
        end
    end
end

fclose(fileID);