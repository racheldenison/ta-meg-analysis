% channel_threshold_tests.m

%% load data
a = load('/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG/R1028_20151216/mat/peakMeansStimAve_30Hz_ebi.mat');
s(2) = load('/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG/R0817_20150504/mat/peakMeansStimAve_30Hz_ebi.mat');
s(3) = load('/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG/R1027_20151216/mat/peakMeansStimAve_30Hz_ebi.mat');

s(1).peakMeansStimAve = a.peakMeansStimAve;
s(1).peakMeansBlankAve = a.peakMeansBlankAve;

nS = numel(s);

%% number of channels passing different thresholds
for iS = 1:nS
    nChBlank5(iS) = nnz(s(iS).peakMeansStimAve > median(s(iS).peakMeansBlankAve) + 5*iqr(s(iS).peakMeansBlankAve));
    nChBlank10(iS) = nnz(s(iS).peakMeansStimAve > median(s(iS).peakMeansBlankAve) + 10*iqr(s(iS).peakMeansBlankAve));
    nChStim3(iS) = nnz(s(iS).peakMeansStimAve > median(s(iS).peakMeansStimAve) + 3*iqr(s(iS).peakMeansStimAve));
    nChStim5(iS) = nnz(s(iS).peakMeansStimAve > median(s(iS).peakMeansStimAve) + 5*iqr(s(iS).peakMeansStimAve));
    nChStim10(iS) = nnz(s(iS).peakMeansStimAve > median(s(iS).peakMeansStimAve) + 10*iqr(s(iS).peakMeansStimAve));
end

% trying other ones
for iS = 1:nS
    nChBlank15(iS) = nnz(s(iS).peakMeansStimAve > median(s(iS).peakMeansBlankAve) + 15*iqr(s(iS).peakMeansBlankAve));
end

nChAll = [nChBlank5; nChBlank10; nChStim3; nChStim5; nChStim10];

%% plot
figure
bar(nChAll)
set(gca,'xticklabel',{'Blank5', 'Blank10', 'Stim3', 'Stim5', 'Stim10'})
legend('good SNR','average SNR','bad SNR')
ylabel('number of channels passing threshold')
