function behav = behavior(behav)

%% key
% blockNames = {'blank','fast-left'}; % fast-left
% attBlockNames = {'no-att','att-right'}; % att-right
% targetBlockNames = {'no-targ','pres-pres','pres-abs','abs-pres','abs-abs'};
% cueBlockNames = {'no-cue','1-1','1-2','2-1','2-2','3-1','3-2'}; % 2-1 = cueT2,postcueT1

%% organize data
nTrials = size(behav.responseData_all,1);

cueCondIdx = strcmp(behav.responseData_labels,'cue condition');
t1CondIdx = strcmp(behav.responseData_labels,'target type T1');
t2CondIdx = strcmp(behav.responseData_labels,'target type T2');
responseIdx = strcmp(behav.responseData_labels,'response');
correctIdx = strcmp(behav.responseData_labels,'correct');
rtIdx = strcmp(behav.responseData_labels,'RT');

cueCond = behav.responseData_all(:,cueCondIdx);
t1Cond = behav.responseData_all(:,t1CondIdx);
t2Cond = behav.responseData_all(:,t2CondIdx);
targetCond = [t1Cond t2Cond];
response = behav.responseData_all(:,responseIdx);
correct = behav.responseData_all(:,correctIdx);
rt = behav.responseData_all(:,rtIdx);

%% cue type (T1, T2)
cuedTarget = nan(nTrials,1);
cuedTarget(cueCond==2 | cueCond==3) = 1; % '1-1','1-2' cue T1
cuedTarget(cueCond==4 | cueCond==5) = 2; % '2-1','2-2' cue T2
cuedTarget(cueCond==6 | cueCond==7) = 3; % '3-1','3-2' cue neutral

%% cue type (valid, invalid)
cueValidity = nan(nTrials,1);
cueValidity(cueCond==2 | cueCond==5) = 1; % '1-1','2-2' valid
cueValidity(cueCond==3 | cueCond==4) = -1; % '1-2','2-1' invalid
cueValidity(cueCond==6 | cueCond==7) = 0; % '3-1','3-2' neutral

%% target type (CCW, CW, absent)
responseTarget = nan(nTrials,1);
responseTarget(cueCond==2 | cueCond==4 | cueCond==6) = 1; % '1-1','2-1','3-1'
responseTarget(cueCond==3 | cueCond==5 | cueCond==7) = 2; % '1-2','2-2','3-2'

targetOrientation = nan(nTrials,1);
for i = 1:nTrials
    if responseTarget(i)~=0 && ~isnan(targetCond(i,1))
        targetOrientation(i) = targetCond(i,responseTarget(i));
    end
end

% sanity check (compare to response and correct)
correctResponse = targetOrientation;
correctResponse(targetOrientation==0) = 3;

%% detection performance
targetPresent = nan(nTrials,1);
targetPresent(targetOrientation==1 | targetOrientation==2) = 1;
targetPresent(targetOrientation==0) = 0;

presentResponse = nan(nTrials,1);
presentResponse(response==1 | response==2) = 1;
presentResponse(response==3) = 0;

detectHit = targetPresent==1 & presentResponse==1;
detectMiss = targetPresent==1 & presentResponse==0;
detectFA = targetPresent==0 & presentResponse==1;
detectCR = targetPresent==0 & presentResponse==0;

detectHMFC = double([detectHit detectMiss detectFA detectCR]);
detectHMFC(isnan(targetPresent),:) = NaN;
detectHMFC(targetPresent==0,1:2) = NaN; % target present trials only for H & M
detectHMFC(targetPresent==1,3:4) = NaN; % target absent trials only for FA & CR

%% discrimination performance
discrimCorrect = targetPresent==1 & correct==1;
discrimIncorrect = targetPresent==1 & presentResponse==1 & correct==-1;

discrimCI = double([discrimCorrect discrimIncorrect]);
discrimCI(isnan(targetPresent),:) = NaN;
discrimCI(targetPresent==0,:) = NaN;

discrimHit = targetOrientation==1 & response==1;
discrimMiss = targetOrientation==1 & response==2;
discrimFA = targetOrientation==2 & response==1;
discrimCR = targetOrientation==2 & response==2;

discrimHMFC = double([discrimHit discrimMiss discrimFA discrimCR]);
discrimHMFC(isnan(targetOrientation),:) = NaN;
discrimHMFC(targetOrientation==2,1:2) = NaN; % target 1 trials only for H & M
discrimHMFC(targetOrientation==1,3:4) = NaN; % target 2 trials only for FA & CR

%% overall accuracy
acc = behav.responseData_all(:,correctIdx);
acc(acc==-1) = 0;

%% exclude missed responses and wrong button presses
wWrongButton = behav.responseData_all(:,responseIdx)==0;
wNotMissed = behav.responseData_all(:,correctIdx)==1 | behav.responseData_all(:,correctIdx)==-1;
w = wWrongButton | ~wNotMissed;
detectHMFC(w,:) = NaN;
discrimCI(w,:) = NaN;
discrimHMFC(w,:) = NaN;
acc(w,:) = NaN;
rt(w,:) = NaN;

%% store
behav.cuedTarget = cuedTarget;
behav.cueValidity = cueValidity;
behav.responseTarget = responseTarget;
behav.targetOrientation = targetOrientation;
behav.targetPresent = targetPresent;
behav.presentResponse = presentResponse;
behav.detectHMFC = detectHMFC;
behav.discrimCI = discrimCI;
behav.discrimHMFC = discrimHMFC;
behav.acc = acc;
behav.rt = rt;

