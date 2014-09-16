function rd_checkTriggers(fileName, trigChans)

% %% setup
% exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
% sessionDir = 'R0890_20140806';
% dataFile = 'R0890_TAPilot_8.06.14.sqd';
% fileName = sprintf('%s/%s/%s', exptDir, sessionDir, dataFile);
% 
% trigChans = 160:166;

%% plot triggers from each channel
figure
hold all

for iTrig = 1:numel(trigChans)
    triggers = all_trigger(fileName, trigChans(iTrig));
    trigTimes = triggers(:,1);

    plot(trigTimes, ones(size(trigTimes))+iTrig-1, '.')
    
    fprintf('Channel %d: %d triggers\n', trigChans(iTrig), numel(trigTimes))
end
ylim([0 10])
xlabel('time')
ylabel('triggers')