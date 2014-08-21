function trigs1 = addBlockCodeToTrigger(trigs0, targetTriggers)
%
% function trigEvents = addBlockCodeToTrigger(trigEvents, targetTriggers)
%
% Let's say you have a trigger that appears within a block (eg. for a
% target event), but you didn't code everything about the block within the
% target trigger. This function goes through your target triggers and looks
% up what block they appear in by finding the previous event that is *not*
% a target -- assuming that is the block trigger.
%
% INPUTS:
% trigs0 is a vector with the original triggers
% targetTriggers are the trigger codes to be updated
%
% OUTPUTS:
% trigs1 is the new trigger vector, which is the original trigger vector +
% 100*block trigger for all target triggers. all non-target triggers remain
% unchanged.
%
% Rachel Denison
% August 2014

% targetTriggers = [16 32];
% trigs0 = trigEvents(:,2);
trigs1 = trigs0;
for iTrig = 1:numel(trigs0)
    if ismember(trigs0(iTrig), targetTriggers)
        blockIdx = find(~ismember(trigs0(1:iTrig),targetTriggers),1,'last');
        trigs1(iTrig) = trigs0(iTrig)+trigs0(blockIdx)*100;
    else
        trigs1(iTrig) = trigs0(iTrig);
    end
end
    