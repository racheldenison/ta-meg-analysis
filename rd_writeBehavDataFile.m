% rd_writeBehavDataFile.m

%% i/o
exptDir = pathToTANoise('Behavior');
analysisDir = sprintf('%s/Group/mat', exptDir);

analysisName = 'gN10_behavior_workspace';

textfileName = 'TANoise_N10_rt'; % dprime, rt

%% load data
D = load(sprintf('%s/%s.mat', analysisDir, analysisName));

%% setup
data = D.groupData.rt; % discrimDprime, rt

nValidity = size(data,1);
nTarget = size(data,2);
nSubSess = size(data,3);

targetNames = {'t1','t2'};
attNames = {'valid','invalid'};

%% write text file
fileID = fopen(sprintf('%s/%s.txt', analysisDir, textfileName),'w');
fprintf(fileID,'%s %s %s %s %s\n','subject','session',...
    'target','att','measure');

    for iSS = 1:nSubSess
        subject = ceil(iSS/2);
        session = 2-mod(iSS,2);
        
        for iT = 1:nTarget
            target = targetNames{iT};
            
            for iV = 1:nValidity
              
                att = attNames{iV};
                
                val = data(iV,iT,iSS);
                
                fprintf(fileID,'%d %d %s %s %1.4f\n', ...
                    subject, session, target, att, val);
            end
        end
    end
    
fclose(fileID);