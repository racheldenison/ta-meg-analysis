function rd_TA2_decode(exptDir, sessionDir)

%% i/o
% exptDir = '/Local/Users/denison/Data/TA2/MEG';
% sessionDir = 'R1187_20181119';

dataFile = dir(sprintf('%s/%s/mat/*condData.mat', exptDir, sessionDir));
dataFileName = sprintf('%s/%s/mat/%s', exptDir, sessionDir, dataFile.name);
figDir = sprintf('%s/%s/figures/ebi_ft', exptDir, sessionDir);

saveFigs = 0;
saveAnalysis = 0;

decodeAll = 1; % decode all trials (1) or by cueing condition (0)
getWeights = 0;
syntheticTrials = 0;

if decodeAll
    analysisFileName = sprintf('%s/%s/mat/classAccT1T2All', exptDir, sessionDir);
else
    analysisFileName = sprintf('%s/%s/mat/classAccT1T2ByCond', exptDir, sessionDir);
end

analStr = 'sp5_nt5';

%% load data
load(dataFileName)

% load data header for plotting topologies
load data/data_hdr.mat

%% setup
t = D.t;
Fs = D.Fs;
eventTimes = D.eventTimes;
data0 = D.condData;

if iscell(data0)
    data = data0;
else
    data = [];
    for iCue = 1:size(data0,4)
        for iT1 = 1:size(data0,5)
            for iT2 = 1:size(data0,6)
                data{iCue,iT1,iT2} = data0(:,:,:,iCue,iT1,iT2);
            end
        end
    end
end

nT = numel(t);
sz = size(data);
nCue = sz(1);
nT1 = sz(2);
nT2 = sz(3);

if decodeAll
    cueNames = {'all trials'};
    figName = {'classAccT1T2All'};
else
    cueNames = {'precue T1','precue T2','neutral'};
    figName = {'classAccT1T2ByCond'};
end

%% remove nan data
% samplingInterval = 1;
% tau = 100;
% filtTau = samplingInterval/tau;

dataRaw = [];
dataFilt = [];
for iCue = 1:nCue
    for iT1 = 1:nT1
        for iT2 = 1:nT2
            nTrials = size(data{iCue,iT1,iT2},3);
            for iTrial = 1:nTrials
                vals = data{iCue,iT1,iT2}(:,:,iTrial);
                idx = isnan(vals(1,:));
                vals(:,idx) = [];
                
                dataRaw{iCue,iT1,iT2}(:,~idx,iTrial) = vals;
                
%                 valsfilt = filter([1-filtTau filtTau-1],[1 filtTau-1], vals);
%                 dataFilt{iCue,iT1,iT2}(:,~idx,iTrial) = valsfilt;
            end
        end
    end
end

%% decoding setup
dataInput = dataRaw; % dataRaw, dataFilt

targetNames = {'T1','T2'};
targetWindows = {[1000 1400],[1300 1700]};
nTarget = numel(targetNames);

nSynTrials = 100; % if constructing synthetic trials
nt = 5; % average this many trials together to improve SNR
sp = 5; % sampling period
kfold = 5;
svmops = sprintf('-s 0 -t 0 -c 1 -v %d -q', kfold);
svmopsNoCV = '-s 0 -t 0 -c 1 -q';

if syntheticTrials
    nReps = 1;
else
    nReps = nt;
end

%% decoding
if decodeAll
    nCue = 1;
    nTarget = 1;
    
%     % grid search
%     nTarget = 1;
%     cParams = 2.^(-1:.5:3); %2.^(-5:2:15);
%     nCParams = numel(cParams);
end

% channels = [15 60 26 14 43 23 26 8 7 1 50 51 2 20 25 13 32 63];
channels = 1:157;

classAccNT = [];
for iRep = 1:nReps
% % grid search
% classAccC = [];
% for iC = 1:nCParams
%     svmops = sprintf('-s 0 -t 0 -c %f -v %d -q', cParams(iC), kfold);
%     disp(svmops)
classAcc = [];
for iT = 1:nTarget
    target = targetNames{iT};
    
    twin = targetWindows{iT};
    times = twin(1):sp:twin(2);
    
    for iCue = 1:nCue
        fprintf('\n%s\n',cueNames{iCue})
        if decodeAll
            switch target
                case 'T1'
                    d1 = dataInput(:,1,:); % T1 vertical
                    d2 = dataInput(:,2,:); % T1 horizontal
                case 'T2'
                    d1 = dataInput(:,:,1); % T2 vertical
                    d2 = dataInput(:,:,2); % T2 horizontal
            end
        else
            switch target
                case 'T1'
                    d1 = dataInput(iCue,1,:); % T1 vertical
                    d2 = dataInput(iCue,2,:); % T1 horizontal
                case 'T2'
                    d1 = dataInput(iCue,:,1); % T2 vertical
                    d2 = dataInput(iCue,:,2); % T2 horizontal
            end
        end
        d1 = d1(:); d2 = d2(:);

        vals1 = []; vals2 = [];
        for i = 1:numel(d1)
            vals1 = cat(3, vals1, d1{i});
            vals2 = cat(3, vals2, d2{i});
        end
        
        % average trials
        if nt > 1
            vals1a = []; vals2a = [];
            n = size(vals1,3);
            if syntheticTrials
                nIdx = nSynTrials*nt;
                trialsIdx = [];
                for i = 1:ceil(nIdx/n)
                    trialsIdx = [trialsIdx randperm(n)];
                end
                startTrials = 1:nt:nIdx;
            else
                trialsIdx = randperm(n);
                startTrials = 1:nt:n;
            end
            for iST = 1:numel(startTrials)
                trIdx = trialsIdx(startTrials(iST):startTrials(iST)+nt-1);
                vals1a(:,:,iST) = mean(vals1(:,:,trIdx),3);
                vals2a(:,:,iST) = mean(vals2(:,:,trIdx),3);
            end
            vals1 = vals1a; vals2 = vals2a;
        end
        
        vals0 = cat(3, vals1, vals2);
        labels0 = [ones(size(vals1,3),1); zeros(size(vals2,3),1)];
        
        %% stratify
        nSamples = numel(labels0);
        foldSize = ceil(nSamples/kfold/2); % 2 classes
        stratIdx = [];
        for iFold = 1:kfold
            idx1 = (1:foldSize) + (iFold-1)*foldSize;
            idx2 = idx1 + nSamples/2;
            stratIdx = [stratIdx idx1 idx2];
        end
        stratIdxS = sort(stratIdx);
        r = stratIdxS(diff(stratIdxS)==0);
        ridx = [];
        for iR = 1:numel(r)
            ridx(iR) = find(stratIdx==r(iR),1,'last');
        end
        stratIdx(ridx) = [];
        if numel(stratIdx)>numel(labels0)
            stratIdx(numel(labels0)+1:end) = [];
        end
        
        vals = vals0(:,channels,stratIdx);
        labels = labels0(stratIdx);
        
        %% classify
        tic
        acc = [];
        for iTime = 1:numel(times)
            fprintf(' ')
            time = times(iTime);
            
            % classification data
            X = squeeze(mean(vals(find(t==time):find(t==time+sp-1),:,:),1))'; % average across time window
            Y = labels;
            
            % remove nan
            idx = isnan(X(:,1));
            X(idx,:) = [];
            Y(idx) = [];
            
            % scale data
            Xs = zscore(X);
%             Xss = Xs./repmat(max(abs(Xs)),size(Xs,1),1); % range [-1,1]
            
            % fit and cross validate classifier
            acc(iTime) = svmtrain(Y, Xs, svmops);
            
%             % example of separate prediction and classification steps
%             model1 = svmtrain(trainlabels1, trainfeatures1, '-s 0 -t 0 -c 1');
%             predlabels = svmpredict(testlabels1, testfeatures1, model1);
%             predacc = mean(predlabels==testlabels1);
            
            % get the svm model, no cv
            if getWeights
                model(iTime) = svmtrain(Y, Xs, svmopsNoCV);
            else
                model = [];
            end
        end
        toc
        
        classAcc(:,iCue,iT) = acc;
        classModel{iCue,iT} = model;
    end
end
% % grid search
% classAccC(:,:,:,iC) = classAcc;
% end

% trial average
classAccNT(:,:,:,iRep) = classAcc;
classModelNT(:,:,iRep) = classModel;
end

%% extract channel weights 
if getWeights
    classWeights = [];
    for iT = 1:nTarget
        twin = targetWindows{iT};
        times = twin(1):sp:twin(2);
        for iCue = 1:nCue
            for iTime = 1:numel(times)
                for iRep = 1:nReps
                    model = classModelNT{iCue,iT,iRep}(iTime);
                    
                    w = model.SVs' * model.sv_coef;
                    b = -model.rho;
                    if (model.Label(1) == -1)
                        w = -w; b = -b;
                    end
                    classWeights(:,iTime,iCue,iT,iRep) = w;
                end
            end
        end
    end
else
    classWeights = [];
end

%% plot
ylims = [30 80];
classTimes = [];

figure
for iT = 1:nTarget
    twin = targetWindows{iT};
    times = twin(1):sp:twin(2);
    xlims = twin;
    
    subplot(nTarget,1,iT)
    hold on
    plot(times, mean(classAccNT(:,:,iT,:),4),'LineWidth',1)
    plot(xlims,[50 50],'k')
    xlim(xlims)
    ylim(ylims)
    xlabel('time (ms)')
    title(sprintf('T%d',iT))
    if iT==1
        legend(cueNames)
    else
        ylabel('classification accuracy (%)')
    end
    
    classTimes(:,iT) = times; % store classification times
end

if saveFigs
    rd_saveAllFigs(gcf, {sprintf('%s_%s',figName{1},analStr)}, 'plot', figDir)
end

%% topo weights movie T1 and T2
if getWeights
    clims = [-3 3];
    
    figure('Position',[250 850 950 450])
    for iTime = 1:size(classTimes,1)
        for iT = 1:nTarget
            vals = squeeze(mean(classWeights(:,iTime,1,iT,:),5))';
            subplot(1,nTarget,iT)
            ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
            set(gca,'CLim',clims)
            colorbar
            title(sprintf('t = %d', classTimes(iTime,iT)))
        end
        pause(0.2)
        %     input('go')
    end
end

%% store results
A.cueNames = cueNames;
A.targetNames = targetNames;
A.targetWindows = targetWindows;
A.decodingOps.channels = channels;
A.decodingOps.nTrialsAveraged = nt;
A.decodingOps.binSize = sp;
A.decodingOps.kfold = kfold;
A.decodingOps.svmops = svmops;
A.classTimes = classTimes;
A.classAcc = classAcc;
A.classModel = classModel;
A.classWeights = classWeights;

%% save analysis
if saveAnalysis
    save(sprintf('%s_%s.mat',analysisFileName,analStr), 'A')
end
