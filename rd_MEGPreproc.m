function preprocFileName = rd_MEGPreproc(filename, figDir)

%% Setup
% desk
% filename = '/Local/Users/denison/Data/TAPilot/MEG/R0817_20140820/R0817_TAPilot_8.20.14.sqd';
% filename = '/Local/Users/denison/Data/TAPilot/MEG/R0817_20140820/preproc/R0817_TAPilot_8.20.14_run01.sqd';
% figDir = '/Local/Users/denison/Data/TAPilot/MEG/R0890_20140806/Runs/figures';

% remember, these channel numbers use zero indexing
megChannels = 0:156;
refChannels = 157:159;
triggerChannels = 160:166;
photodiodeChannel = 191;

% badChannels = [];

% preproc options
Fl = 60; % line noise frequency
environmentalDenoise = 1;
applyLineNoiseFilter = 0;
removeBadChannels = 1;
TSPCA = 0;
interpolate = 0;

plotFigs = 1;
saveFigs = 1;

analStr = [];

%% Get the MEG data
% data is time x channels
[data, info] = sqdread(filename);

Fs = info.SampleRate;
t = 0:1/Fs:size(data,1)/Fs-1/Fs;

%% Set aside data from special channels
% add 1 to adjust for zero-indexing
refData = data(:,refChannels+1);
trigData = data(:,triggerChannels+1);
pdData = data(:,photodiodeChannel+1);

%% Look at special channels
if plotFigs
    figure
    subplot(3,1,1)
    plot(t,refData)
    title(['reference channels' num2str(refChannels)])
    subplot(3,1,2)
    plot(t,trigData)
    legend(num2str(triggerChannels'))
    title(['trigger channels' num2str(triggerChannels)])
    subplot(3,1,3)
    plot(t,pdData)
    title(['photodiode channel' num2str(photodiodeChannel)])
    xlabel('time (s)')
end

%% Denoise using reference channels
if environmentalDenoise
    % See also LSdenoise.m
    analStr = [analStr 'e'];
    
    % convert to time x trials x channels
    data = permute(data,[1 3 2]);
    data = meg_environmental_denoising(data, refChannels+1, megChannels+1, plotFigs);
    data = permute(data,[1 3 2]); % convert back
end

%% Line noise filter
if applyLineNoiseFilter
    analStr = [analStr 'l'];
    
    % data should be channels x time
    data = ft_preproc_dftfilter(data', Fs, Fl);
    data = data';
end

%% Find bad channels
if removeBadChannels
    analStr = [analStr 'b'];
    
    % high or low variance across the entire time series
    outlierSDChannels = meg_find_bad_channels(permute(data(:,megChannels+1),[1 3 2]));
    
    % dead or saturating channels for all or portions of the time series
    deadChannels = checkForDeadChannels(filename)+1;
    
    % aggregate the bad channels
    badChannels = unique([outlierSDChannels' deadChannels]);
    nBad = numel(badChannels);
    
    % plot the time series for the bad channels
    if plotFigs
        figure
        for iBad = 1:nBad
            subplot(nBad,1,iBad)
            chan = badChannels(iBad);
            plot(t, data(:,chan))
            title(sprintf('channel %d', chan))
        end
        xlabel('time (s)')
    end
    
    % zero out bad channels
    data(:,badChannels) = 0;
end

%% Save data preprocessed up to this point
preFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);

if exist(preFile,'file')
    error('preFile already exists ... will not overwrite. exiting.')
else
    sqdwrite(filename, preFile, 'data', data);
end

%% Time-shift PCA for environmental denoising
% http://www.isr.umd.edu/Labs/CSSL/simonlab/Denoising.html
% http://lumiere.ens.fr/Audition/adc/meg/
% "Noise fields measured by reference magnetometers are optimally filtered 
% and subtracted from brain channels. The filters (one per reference/brain 
% sensor pair) are obtained by delaying the reference signals, 
% orthogonalizing them to obtain a basis, projecting the brain sensors onto 
% the noise-derived basis, and removing the projections to obtain clean 
% data." (DOI: 10.1016/j.jneumeth.2007.06.003)
if TSPCA
    sizeOfBlocks = 20000;
    shifts = -100:100;
    sourceFile = preFile;
    
    analStr = [analStr 't'];
    tspcaFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    
    % run sqd denoise
    % this writes the tspca sqd file
    if exist(tspcaFile,'file')
        error('tspcaFile already exists ... will not overwrite. note that continuing the script will delete this file.')
    else
        fprintf('Running sqdDenoise\n');
        sqdDenoise(sizeOfBlocks, shifts, 0, sourceFile, badChannels-1, 'no', ...
            168, 'no', tspcaFile);
    end
    
    % plot comparison with basic environmental denoising
    if plotFigs
        datats = sqdread(tspcaFile);
        figure
        sampleChannels = [1 14 badChannels];
        for iCh = 1:numel(sampleChannels)
            chan = sampleChannels(iCh);
            subplot(numel(sampleChannels),1,iCh)
            hold on
            plot(t, data(:,chan))
            plot(t, datats(:,chan), 'r')
            xlabel('time (s)')
            title(sprintf('channel %d', chan))
        end
        legend('LS environmental denoise','+time-shift PCA environmental denoise')
    end
    
    % replace data with datats
    data = datats;
    clear datats;
    
    % remove just-created sourceFile
    delete(sourceFile);
end

%% Interpolate to replace bad channels
if interpolate
    % create dummy ft_data from the original data and update it
    dataset = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    analStr = [analStr 'i'];
    
    % read data into ft structure in continuous mode by initializing cfg with
    % only the dataset
    % ft data is channels x time
    cfg = [];
    cfg.dataset = dataset;
    ft_data = ft_preprocessing(cfg); % this preserves the data, but scales it by like 10^-13
    ft_data.trial{1} = data'; % so just replace data with the original data
    
    % % interpolate bad channels interpolation
    % % if using 'nearest' or 'average' interpolation, get the neighbours
    % cfg = [];
    % cfg.method = 'distance';
    % [neighbours, cfg] = ft_prepare_neighbours(cfg, ft_data);
    
    cfg = [];
    cfg.badchannel = ft_channelselection(badChannels, ft_data.label);
    cfg.method = 'spline';
    % cfg.neighbours = neighbours; % only needed for nearest and average, not for spline
    ft_data = ft_channelrepair(cfg, ft_data);
    
    % compare before and after interpolation
    if plotFigs
        figure
        sampleChannels = badChannels;
        for iCh = 1:numel(sampleChannels)
            chan = sampleChannels(iCh);
            subplot(numel(sampleChannels),1,iCh)
            hold on
            plot(t, data(:,chan))
            plot(t, ft_data.trial{1}(chan,:), 'r')
            xlabel('time (s)')
            title(sprintf('channel %d', chan))
        end
        legend('before interpolation','after interpolation')
    end
    
    % save the interpolated data
    data(:,1:numel(megChannels)) = ft_data.trial{1}';
    
    interpFile = sprintf('%s_%s.sqd', filename(1:end-4), analStr);
    sqdwrite(filename, interpFile, 'data', data);
    
    % if we've run the interpolation, then remove just-created tspcaFile.
    % otherwise, delete the preFile.
    if exist(tspcaFile,'file')
        delete(tspcaFile);
    else
        delete(preFile);
    end
end

%% finally, check the triggers
rd_checkTriggers(filename);

%% save figs
if saveFigs
    runStr = rd_getTag(filename);
    figSubDir = sprintf('%s/%s/%s', figDir, analStr, runStr);
    if ~exist(figSubDir,'dir')
        mkdir(figSubDir)
    end
    f = sort(findobj('Type','figure'));
    for iF = 1:numel(f)
        figNames{iF} = num2str(f(iF));
    end
    rd_saveAllFigs(f,figNames,[],figSubDir);
    close all
end

%% return preproc file name
preprocFileName = sprintf('%s_%s.sqd', filename(1:end-4), analStr);




