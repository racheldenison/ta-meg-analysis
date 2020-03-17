function [ft_cleandata, WorkFlow, pcStr, ft_PCA] = meg_pca_ica(dataset, badChannels, trialDef, rejectPC, rejectIC)
%
% function [ft_cleandata, WorkFlow] = meg_pca_ica(dataset, badChannels, trialDef)
%
% Performs PCA then ICA, allowing the manual rejection of components at
% each stage. PCA first sorts components according to variance, and then 
% ICA looks only at the first 32 components.
%
% Returns a fieldtrip structure with the clean data (reconstructed after
% rejecting components) stored in the 'trial' field. Also returns a
% structure WorkFlow, which has a record of the steps performed.
%
% Rachel Denison
% Dec 2014
% Based largely on code from Adeen Flinker, MEG_AnalysisFlow_SingleTone.m
%
%% sample inputs
% dataset = 'R0890_Singletone_8.06.14.sqd';
% badChannels = [];
% 
% trialDef.trialFunHandle = @mytrialfun_all;
% trialDef.prestim = 0.5;
% trialDef.poststim = 1;
% trialDef.trig = 161:167;
% trialDef.nTrigsExpected = 100;

%% add paths
addpath('/e/1.3/p1/denison/Software/fieldtrip-20170510/external/eeglab')

%% store the inputs
WorkFlow.dataset = dataset;
WorkFlow.badChannels = badChannels;
WorkFlow.trialDef = trialDef;

%% put data into fieldtrip format
% data channels
cfg = struct('dataset',dataset,'channel',[1:157],'demean','no');
ft_data = ft_preprocessing(cfg);
data = ft_data.trial{1};

% remove bad-channels
data(badChannels,:) = 0;

% extract trials
cfg.trialdef.prestim    = trialDef.prestim;
cfg.trialdef.poststim   = trialDef.poststim;
cfg.trialdef.trig       = trialDef.trig;
threshold               = trialDef.threshold;
nTrigsExpected          = trialDef.nTrigsExpected;

[trl,Events] = trialDef.trialFunHandle(cfg, threshold, nTrigsExpected);

% data inclusion
WorkFlow.data_continuous_block = [1 length(data)];
WorkFlow.Events = Events;

pcStr = []; 

%% run pca
% Remove mean and transpose to channels x samples
DataDemean = data - repmat(mean(data,2),1,size(data,2));
[EigenVectors,EigenValues]=pcsquash(DataDemean);            % pca eigenvectors from eeglab
WorkFlow.PCAmixing = EigenVectors';

ft_data.trial{1} = DataDemean; %is the whole experiment. 
ft_PCA = ft_componentanalysis(struct('demean','no','unmixing',EigenVectors','topolabel',{ft_data.label}),ft_data);

% view PCA
layout = ft_prepare_layout(ft_data.cfg,ft_data);
WorkFlow.layout = layout;

fs = 1000;
NSec = 4;
figPos = [1 5 800 1364];
figure
Aft_plot_component_rd(ft_PCA,1:5,layout,trl,WorkFlow.data_continuous_block(1)-1,fs,NSec,figPos);

saveas(gcf,'output','jpg');
WorkFlow.PCA_screenshot = imread('output.jpg');
%ft_databrowser(struct('viewmode','component','layout',layout), ft_PCA);

%  reject components (or  not)
% reject_PCA_comps = input('\nInput PCA components to reject (e.g. [1 4] or []): ');
if rejectPC
    reject_PCA_comps = 1;
else
    reject_PCA_comps = 0; 
end
fprintf('\nAutomatically rejecting 1st PC\n')
if ~isempty(reject_PCA_comps)
    PCA_postreject = ft_rejectcomponent(struct('component',reject_PCA_comps,'demean','no'),ft_PCA);
    %retun PCA after rejection
    DataDemean = PCA_postreject.trial{1} - repmat(mean(PCA_postreject.trial{1},2),1,size(PCA_postreject.trial{1},2));
    [EigenVectors,EigenValues]=pcsquash(DataDemean);            % pca eigenvectors from eeglab
    
    ft_data.trial{1} = DataDemean;
    ft_PCA = ft_componentanalysis(struct('demean','no','unmixing',EigenVectors','topolabel',{ft_data.label}),ft_data);
    pcStr = [pcStr 'p']; 
end

WorkFlow.PCA_rejected_components = reject_PCA_comps;
WorkFlow.postrejection_PCAmixing = EigenVectors';

%% run ica
ncomps = 32;  % number of subspace PCA components
[weights, sphere] = runica(ft_PCA.trial{1}(1:ncomps,:),'lrate',0.001);

dummy = ft_PCA; dummy.trial{1} = ft_PCA.trial{1}(1:ncomps,:); dummy.topo = dummy.topo(1:ncomps,1:ncomps);
dummy.label = dummy.label(1:ncomps); dummy.topolabel = dummy.topolabel;

ft_ICA = ft_componentanalysis(struct('demean','no','unmixing',weights*sphere,'topolabel',{dummy.label}),dummy);

% change topography data to match ica + pca transformation
ft_ICA.topolabel_orig = ft_ICA.topolabel;
ft_ICA.topolabel = ft_PCA.topolabel;
ft_ICA.topo = pinv(weights*sphere*EigenVectors(:,1:ncomps)'*eye(length(EigenVectors)));

WorkFlow.ICA_PCAncomps = ncomps;
WorkFlow.ICAmixing = weights*sphere;
WorkFlow.ICA_topo = ft_ICA.topo;

%% view ica
% event related averaged ICA and frequency spectrum

% ft_databrowser(struct('viewmode','component','layout',layout), ft_ICA);
% close all;
figPos = [1 5 800 1364];
Nsec = 4;
ER_shift = WorkFlow.data_continuous_block(1)-1;
compSets = {1:6,7:12,13:18,19:24,25:30,31:32};
latestFig = gcf;
for i = numel(compSets):-1:1
    figure(i+latestFig.Number)
    Aft_plot_component_rd(ft_ICA,compSets{i},layout,trl,ER_shift,fs,Nsec,figPos);
end

%%  ica reject components (or  not)
% select and visualize components to reject
if rejectIC
    reject_ICA_comps = input('\nInput ICA components to reject (e.g. [1 4] or []): ');
else 
    reject_ICA_comps = [];
end
fprintf('\nAutomatically not rejecting any ICs\n')

if ~isempty(reject_ICA_comps)
    figure
    Aft_plot_component_rd(ft_ICA,reject_ICA_comps,layout,trl,ER_shift,1000,Nsec,figPos);
    saveas(gcf,'output','jpg');
    WorkFlow.ICA_rejected_screenshot = imread('output.jpg');
    pcStr = [pcStr 'c'];
end

WorkFlow.ICA_rejected_components = reject_ICA_comps;

% backprojection
activations = ft_ICA.trial{1};
activations(reject_ICA_comps,:) = 0;  % remove components
ICA_postreject = inv(weights*sphere)*activations;

ft_PCA_ICA = ft_PCA;
ft_PCA_ICA.trial{1}(1:ncomps,:) = ICA_postreject;
ft_cleandata = ft_rejectcomponent(struct('component',[],'demean','no'),ft_PCA_ICA);

%% view cleandata (blue) and original data (red)
% windowSize = [1 5 2560 1392];
% eegplot(ft_cleandata.trial{1}./1e-13,'srate',ft_cleandata.fsample,'winlength',5,'dispchans',50,'position',windowSize,'data2',data./1e-13);


