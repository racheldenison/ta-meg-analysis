% fixTrialsRejected.m

% load trials_rejected
trialsRejected0 = trials_rejected;

ftTrials = zeros(1,480);
ftTrials(trialsRejected0) = 1;

trials = ones(1,516);
b = 1:21:43;
blankTrials = [];
for i = 1:12
    blankTrials = [blankTrials (i-1)*43+b];
end
trials(blankTrials) = nan;

trials(trials==1) = ftTrials;

trials_rejected = find(trials==1)';

save('~/Desktop/trials_rejected.mat','trials_rejected')