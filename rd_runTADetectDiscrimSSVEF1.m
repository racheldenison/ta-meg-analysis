% rd_runTADetectDiscrimSSVEF1.m

%% setup
subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
    'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
    'R1019_20151118','R1021_20151120','R1026_20151211', ...
    'R0852_20151211','R1027_20151216','R1028_20151216',...
    'R1029_20151222'}; % N=16
nSubjects = numel(subjects);

%% run analysis
for iSubject = 1:nSubjects
    % get fileBase
    sessionDir = subjects{iSubject};
    fprintf('%s\n',sessionDir)
    rd_TADetectDiscrimSSVEF1(sessionDir);
    close all;
end
fprintf('done.\n')
    
    
    
    
    
    