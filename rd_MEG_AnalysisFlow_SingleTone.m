% MEG_ANALYSISFLOW_SINGLETONE.m
%
% dependencies
%      fieldtrip  http://fieldtrip.fcdonders.nl/download
%      sqddeniose http://www.isr.umd.edu/Labs/CSSL/simonlab/Resources.html
%      sqdproject http://www.isr.umd.edu/Labs/CSSL/simonlab/Resources.html
%      eeglab for eegplot and dependencies 
%% Define directory/dataset
clear all
dr = dir('*ing*one*.sqd')
dataset = dr(1).name;
WorkFlow.dataset = dataset;

%% load/view raw data
data = sqdread(dataset);    % available at http://www.isr.umd.edu/Labs/CSSL/simonlab/Resources.html (sqdproject)
data  = data(:,1:157)';
srate = 1000;
windowSize = [1           5        2560        1392];
%data_view =   ft_resampledata(struct('resamplefs', 40,'detrend','no','demean','no'),data);

%% view data
eegplot(data,'srate',srate,'winlength',20,'dispchans',80,'position',windowSize);
fprintf('dead channels:\n');
deadchannels = checkForDeadChannels(dataset)+1  % available at http://www.isr.umd.edu/Labs/CSSL/simonlab/Resources.html (Denoising)
%% mark bad channels and denoise raw data
clear data;
badchannels = unique([57 deadchannels]);
tspca = 1;

WorkFlow.denoise_command = @LSdenoise;
% WorkFlow.denoise_params = {'~/Desktop/Adeen_test_room_rest.sqd', dataset,[],[],badchannels,[],1};
WorkFlow.denoise_params = {WorkFlow.dataset, dataset,[],[],badchannels,[],1}; % always do tspca, will save both files

WorkFlow.denoise_command(WorkFlow.denoise_params{:});

if tspca
    dataset = [WorkFlow.dataset(1:end-4) '-LSdenoised_NR.sqd']; % update dataset name to denoised data
else
    dataset = [WorkFlow.dataset(1:end-4) '-LSdenoised.sqd']; % update dataset name to denoised data
end
WorkFlow.deadchannels = deadchannels;
WorkFlow.badchannels = badchannels;
WorkFlow.dataset_denoised = dataset;

%% extract desnoised data
% analog channel 191
AnalogChannel = sqdread(dataset,'channels',190);
% data channels
cfg = struct('dataset',dataset,'channel',[1:157],'demean','no');
ft_data = ft_preprocessing(cfg);
data = ft_data.trial{1};

% remove bad-channels
data(badchannels,:) = 0;

% extract triggers and fieldtrip stamps

% singletone sort events - write function per task
% [trl,Events] = SingleTonetrialsort(struct('dataset',dataset,'trialdef',struct('prestim',0.5,'poststim',1,'trig',161:167)));
[trl,Events] = mytrialfun_all(struct('dataset',dataset,'trialdef',struct('prestim',0.5,'poststim',1,'trig',161:167)), 2.5, 100);

% get rid of trials that are too early (RD)
% exclIdx = find(trl(:,1)<0);
% trl(exclIdx,:) = [];
% Events(exclIdx) = [];

%data inclusion
% st_block = 1; 
% en_block = trl(end,2)+min(2000,length(data)-trl(end,2));
% 
% data = data(:,st_block:en_block);
% WorkFlow.TrialArray = trl;
% WorkFlow.data_continuous_block = [st_block en_block];
WorkFlow.data_continuous_block = [1 length(data)];
WorkFlow.Events = Events;
%% run pca
% Remove mean and transpose to channels x samples
DataDemean = data - repmat(mean(data,2),1,size(data,2));
[EigenVectors,EigenValues]=pcsquash(DataDemean);            % pca eigenvectors from eeglab

ft_data.trial{1} = DataDemean;
ft_PCA = ft_componentanalysis(struct('demean','no','unmixing',EigenVectors','topolabel',{ft_data.label}),ft_data);
WorkFlow.PCAmixing = EigenVectors';
% view PCA
layout = ft_prepare_layout(ft_data.cfg,ft_data);
WorkFlow.layout = layout;

figure(1);Aft_plot_component_rd(ft_PCA,1:5,layout,trl,WorkFlow.data_continuous_block(1)-1,1000,3,[1 5 800 1364]);

saveas(gcf,'output','jpg');
WorkFlow.PCA_screenshot = imread('output.jpg');
%ft_databrowser(struct('viewmode','component','layout',layout), ft_PCA);

%  reject components (or  not)
reject_PCA_comps = [3];
if ~isempty(reject_PCA_comps)
    PCA_postreject = ft_rejectcomponent(struct('component',reject_PCA_comps,'demean','no'),ft_PCA);
    %retun PCA after rejection
    DataDemean = PCA_postreject.trial{1} - repmat(mean(PCA_postreject.trial{1},2),1,size(PCA_postreject.trial{1},2));
    [EigenVectors,EigenValues]=pcsquash(DataDemean);            % pca eigenvectors from eeglab
    
    ft_data.trial{1} = DataDemean;
    ft_PCA = ft_componentanalysis(struct('demean','no','unmixing',EigenVectors','topolabel',{ft_data.label}),ft_data);
end

WorkFlow.PCA_rejected_components = reject_PCA_comps;
WorkFlow.postrejection_PCAmixing = EigenVectors';

%% run ica
ncomps = 32;  % number of subspace PCA components
[weights, sphere]=runica(ft_PCA.trial{1}(1:ncomps,:),'lrate',0.001);

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
%  event related averaged ICA and frequency spectrum

%     ft_databrowser(struct('viewmode','component','layout',layout), ft_ICA);
close all;
Pos1 =[ 1         5        800        1364];
Pos2 =[801      5       800        1364];
Pos3 = [1602   5       800        1364];
Nsec = 4;
ER_shift = WorkFlow.data_continuous_block(1)-1;
% figure(1);Aft_plot_component_rd(ft_ICA,[1:10],layout,trl,ER_shift,1000,Nsec,Pos1);
% figure(2);Aft_plot_component_rd(ft_ICA,[11:20],layout,trl,ER_shift,1000,Nsec,Pos2);
% figure(3);Aft_plot_component_rd(ft_ICA,[21:32],layout,trl,ER_shift,1000,Nsec,Pos3);
compSets = {1:6,7:12,13:18,19:24,25:30,31:32};
for i = numel(compSets):-1:1
    figure(i)
    Aft_plot_component_rd(ft_ICA,compSets{i},layout,trl,ER_shift,1000,Nsec,Pos1);
end

%%
reject_ICA_comps = [4 5];
figure(3);Aft_plot_component_rd(ft_ICA,reject_ICA_comps,layout,trl,ER_shift,1000,Nsec,Pos1);
saveas(gcf,'output','jpg');

WorkFlow.ICA_rejected_components = reject_ICA_comps;
WorkFlow.ICA_rejected_screenshot = imread('output.jpg');
%%  ica reject components (or  not)
% backprojection
activations = ft_ICA.trial{1};
activations(reject_ICA_comps,:) = 0;  % remove components
ICA_postreject = inv(weights*sphere)*activations;

ft_PCA_ICA =ft_PCA;
ft_PCA_ICA.trial{1}(1:ncomps,:) = ICA_postreject;
cleandata = ft_rejectcomponent(struct('component',[],'demean','no'),ft_PCA_ICA);

WorkFlow.ICA_rejected_components = reject_ICA_comps;

% view cleandata (blue) and original data (red)
eegplot(cleandata.trial{1}./1e-13,'srate',srate,'winlength',5,'dispchans',50,'position',windowSize,'data2',ft_data.trial{1}./1e-13);
%%
%%   save data with shift and update trial info
%%

% update trl and Event trigger info relative to shift and new data set

shift = WorkFlow.data_continuous_block(1)-1;

%singletone
for i = 1:length(trl)
    Events(i).trigger = trl(i,1)-trl(i,3)-shift;
    trl(i,[1 2]) = trl(i,[1 2])-shift;
end

WorkFlow.Events = Events;
WorkFlow.trl = trl;

data = cleandata.trial{1};

srate = 1000;
dataset  = WorkFlow.dataset; tmp =strfind(WorkFlow.dataset,'_');
dataset = dataset(1:tmp(2));
save([dataset 'preproc'], 'data', 'WorkFlow', 'Events', 'trl', 'srate');


%% downsample, run TF analysis and save data
clear all;
dr= dir('*ing*preproc*.mat');
load(dr(1).name);

new_srate = 100;
[p,q] = rat(new_srate/srate);
data_dn=resample(data',p,q)';

% update trl and Event trigger info
%     for i = 1:length(Events)
%             flds = fieldnames(Events)';
%             for f = flds([1:3 5 7 8])
%                 eval(['Events(i).' f{1} ' = round(Events(i).' f{1} '*p/q);']);
%             end
%             trl(i,[1 2 6]) = round(trl(i,[1 2 6])*p/q);
%     end
%
%  singletone
for i = 1:length(Events)
    flds = fieldnames(Events)';
    for f = flds([1])
        eval(['Events(i).' f{1} ' = round(Events(i).' f{1} '*p/q);']);
    end
    trl(i,[1 2 3]) = round(trl(i,[1 2 3])*p/q);
end


% save data
srate = new_srate;
data = data_dn;


% prepare data
[data_tf, freqs] = Aft_createTF(data,srate,50);

%singletone
evs = [Events.trigger];

close all; Aft_ReviewTrials_rd(data,evs,data_tf,freqs,[],[],[],[],[],WorkFlow.layout);
%% update events
for i = 1:length(Events)
    if ismember(i,WorkFlow.ManualTrials.BadTrials)
        Events(i).isbad =1;
    else
        Events(i).isbad =0;
    end
end

%
evs = [Events.trigger];
isbad = [Events.isbad];
evs = evs(~isbad);
[ER_avg,tm,tm_ms] = Aft_ERaverage(data,evs,-20,100,100,[-20 0]);
[TF_avg,tm,tm_ms] = Aft_ERaverage(data_tf,evs,-20,100,100,[-20 0]);


dataset  = WorkFlow.dataset; tmp =strfind(WorkFlow.dataset,'_');
dataset = dataset(1:tmp(2));
badchannels = union(WorkFlow.badchannels,WorkFlow.ManualTrials.BadChannels);

%% optional reapir bad channels

% create dummy ft_data from the original data and update it
ft_data = ft_preprocessing(struct('dataset',WorkFlow.dataset,'channel',[1:157],'demean','no'));
ft_data.trial{1}=data;
ft_data.time{1}=ft_data.time{1}(1:length(data));
ft_data.sampleinfo = [1 length(data)];
ft_data.cfg.trl = [1 length(data) 0];
ft_data.fsample = new_srate;
% interpolate bad channels using a spline interpolation
interp = ft_channelrepair(struct('method','spline','badchannel',{ft_channelselection(badchannels,ft_data.label)}),ft_data);
WorkFlow.badchannel_interp = 'spline';
data = interp.trial{1};

%% save data
save([dataset 'data'], 'data' ,'srate', 'Events','trl', '*_avg*', 'freqs', 'WorkFlow' ,'badchannels','tm_ms');
save([dataset 'Events'],'Events' ,'trl');
% copyfile([dataset 'data.mat'],['../SingleTone/' dataset 'data.mat']);
