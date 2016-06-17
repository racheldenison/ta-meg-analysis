% rd_plotChannelSelection

%% Setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
analStr = 'ebi'; % '', 'ebi', etc.
ssvefFreq = 40;

subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
    'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
    'R1019_20151118','R1021_20151120','R1026_20151211', ...
    'R0852_20151211','R1027_20151216','R1028_20151216',...
    'R1029_20151222'}; % N=16

nSubjects = numel(subjects);

%% Get data
for iSubject = 1:nSubjects
    sessionDir = subjects{iSubject};
    
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    
    d = load(sprintf('%s/channels_%dHz_ebi.mat', matDir, ssvefFreq));
    
    channelsRanked(:,iSubject) = d.channelsRanked;
    channelsRankedAmps(:,iSubject) = d.channelsRankedAmps;
end
    
%% Plot figs
colors = varycolor(nSubjects);

figure
set(gca,'ColorOrder',colors)
hold all
plot(channelsRankedAmps, 'LineWidth', 2)
legend(subjects)
xlabel('channel (rank-ordered)')
ylabel(sprintf('%d hz amplitude, stim average', ssvefFreq))
% xlim([0 10])

figure
plot(channelsRanked)
xlim([0 20])
    
    
    
    