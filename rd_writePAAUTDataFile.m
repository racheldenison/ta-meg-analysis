% rd_writePAAUTDataFile.m

%% i/o
exptDir = pathToTANoise('MEG');
analysisDir = sprintf('%s/Group/mat', exptDir);

analysisName = 'gN10_itpcVH_20Hz_workspace';

textfileName = 'TANoise_N10_itpc_20Hz';

%% load data
D = load(sprintf('%s/%s.mat', analysisDir, analysisName));

%% setup
data = D.groupData.PAAUT;
t = D.tt;

nTime = size(data,1);
nPAAU = size(data,2);
nTarget = size(data,3);
nSubSess = size(data,4);

targetNames = {'t1','t2'};
vhNames = {'v','h'};
attNames = {'att','unatt'};

%% write text file
fileID = fopen(sprintf('%s/%s.txt', analysisDir, textfileName),'w');
fprintf(fileID,'%s %s %s %s %s %s %s\n','time','subject','session',...
    'target','ori','att','measure');

for iTime = 1:10:nTime
    time = t(iTime);
    for iSS = 1:nSubSess
        subject = ceil(iSS/2);
        session = 2-mod(iSS,2);
        
        for iT = 1:nTarget
            target = targetNames{iT};
            
            for iPAAU = 1:nPAAU
                ori = vhNames{ceil(iPAAU/2)};
                att = attNames{2-mod(iPAAU,2)};
                
                val = data(iTime,iPAAU,iT,iSS);
                
                fprintf(fileID,'%d %d %d %s %s %s %1.4f\n', ...
                    time, subject, session, target, ori, att, val);
            end
        end
    end
end
    
fclose(fileID);