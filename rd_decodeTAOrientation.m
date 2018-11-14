% rd_decodeTAOrientation.m

%% load data
load('/Local/Users/denison/Data/TANoise/MEG/R1103_20180213/mat/analysis_singleTrials_R1103_TANoise_2.13.18_ebi_ft_topChannels5_allTrials_20Hz.mat')

%% setup
t = A.t;
trigNames = A.trigNames;
measure = 'trigMean'; %'wITPC';

twins = {[1000 1400],[1300 1700]};
attConds = {[1 3 5 7],[2 4 6 8]}; % attT1, attT2
oriConds{1} = {[1 2 5 6],[3 4 7 8]}; % T1 V, H
oriConds{2} = {[1 2 3 4],[5 6 7 8]}; % T2 V, H

%% get data
data = A.(measure);

%% preprocessing
switch measure
    case 'trigMean'
        % high-pass filter single trial data to detrend
        samplingInterval = 1;
        tau = 100;
        filtTau = samplingInterval/tau;
        
        nCh = size(data,2);
        nConds = size(data,4);
        dataFilt = nan(size(data));
        
        for iCh = 1:nCh
            for iCond = 1:nConds
                y = squeeze(data(:,iCh,:,iCond));
                idx = isnan(y(1,:));
                y(:,idx) = [];
                yfilt = filter([1-filtTau filtTau-1],[1 filtTau-1], y);
                nGoodTrials = nnz(~idx);
                % yMean(:,iCond) = mean(yfilt,2);
                dataFilt(:,iCh,1:nGoodTrials,iCond) = yfilt;
            end
        end
        
        data0 = data;
        data = dataFilt;
end

%% simple correlation (on wITPC)
for iT = 1:2
    oriCond = oriConds{iT};
    twin = twins{iT};
    tidx = find(t==twin(1)):find(t==twin(2));
    vals = data(tidx,:);
    
    for iAtt = 1:2
        attCond = attConds{iAtt};
        
        classConds{1} = intersect(oriCond{1}, attCond); % vertical
        classConds{2} = intersect(oriCond{2}, attCond); % horizontal
        
        % is correlation higher within class than between classes?
        for iC = 1:numel(classConds)
            rW(iC,iAtt,iT) = corr(vals(:,classConds{iC}(1)),vals(:,classConds{iC}(2)));
        end
        for i = 1:numel(classConds{1})
            rB(i,iAtt,iT) = corr(vals(:,classConds{1}(i)),vals(:,classConds{2}(i)));
        end
    end
end

%% correlation summary
rWMean = squeeze(mean(rW,1));
rBMean = squeeze(mean(rB,1));

figure
hold on
plot(rWMean')
plot(rBMean','--')

%% assemble classification data (on trigMean)
classData = [];
classLabels = [];
for iT = 1:2
    oriCond = oriConds{iT};
    twin = twins{iT};
    tidx = find(t==twin(1)):find(t==twin(2));
    
    for iAtt = 1:2
        attCond = attConds{iAtt};
        
        classConds{1} = intersect(oriCond{1}, attCond); % vertical
        classConds{2} = intersect(oriCond{2}, attCond); % horizontal
    
        vals1 = []; vals2 = [];
        for i = 1:numel(classConds{1})
            vals1 = cat(3, vals1, data(tidx,:,:,classConds{1}(i))); % vertical
            vals2 = cat(3, vals2, data(tidx,:,:,classConds{2}(i))); % horizontal
        end
        
        vals = cat(3, vals1, vals2);
        labels = [ones(size(vals1,3),1); zeros(size(vals2,3),1)];
        
        classData{iAtt,iT} = vals; 
        classLabels{iAtt,iT} = labels;
    end
    
end
    
%% classifier
times = 195:205;
% times = 150:10:350;
% times = 1:size(classData{1,1},1);
nTime = numel(times);


classAcc = [];
iT = 1; iAtt = 1;
fprintf('\nClassifying ... start time: %s\n', datestr(now, 13))
for iTime = 1:nTime
    time = times(iTime);
    if mod(iTime,10)==0
        fprintf('%d, %s\n',iTime, datestr(now, 13))
    else
        fprintf('.')
    end
    X = squeeze(classData{iAtt,iT}(time,:,:))';
    Y = classLabels{iAtt,iT};
    
    % remove nan
    idx = isnan(X(:,1));
    X(idx,:) = [];
    Y(idx) = [];
    
    % fit and cross validate classifier
    svm = fitcsvm(X,Y);
    cv = crossval(svm);
    classAcc(iTime) = 1 - kfoldLoss(cv);
end
    