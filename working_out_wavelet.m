powspctrm     = zeros(nchan,nfoi,ntoi,cfg.precision);

acttboi  = squeeze(~isnan(spectrum(1,1,foiind(ifoi),:)));
nacttboi = sum(acttboi);
if ~hastime
    acttboi  = 1;
    nacttboi = 1;
elseif sum(acttboi)==0
    %nacttboi = 1;
end
acttap = logical([ones(ntaper(ifoi),1);zeros(size(spectrum,1)-ntaper(ifoi),1)]);

powdum = abs(spectrum(acttap,:,foiind(ifoi),acttboi)) .^2;

powspctrm(:,ifoi,acttboi) = powspctrm(:,ifoi,acttboi) + (reshape(mean(powdum,1),[nchan 1 nacttboi]) ./ ntrials);

powspctrm = powspctrm.*ntrials;
powspctrm = powspctrm./trlcnt(ones(size(powspctrm,1),1),:,:);


%%% edited versin
nchan = 1;
nfoi = numel(freqoi);
ntoi = numel(timeoi);
ntaper = 1; ifoi = 1;
ntrials = 1;

powspctrm     = zeros(nchan,nfoi,ntoi);

acttboi  = squeeze(~isnan(spectrum(1,freqIdx,:)));
nacttboi = sum(acttboi);

acttap = logical([ones(ntaper(ifoi),1);zeros(size(spectrum,1)-ntaper(ifoi),1)]);

powdum = abs(spectrum(:,freqIdx,acttboi)) .^2;

powspctrm(:,freqIdx,acttboi) = powspctrm(:,freqIdx,acttboi) + (reshape(mean(powdum,1),[nchan 1 nacttboi]) ./ ntrials);

powspctrm = powspctrm.*ntrials;
powspctrm = powspctrm./trlcnt(ones(size(powspctrm,1),1),:,:);

