function [targetTimes nTargs] = createFixTargetSequence(runDur,frameDur,minDelay,maxDelay)
% runDur = duration of run in seconds
% frameDur = durate of each stimulus frame
% minDelay = minimum target-to-target SOA
% maxDelay = maximum target-to-target SOA

delays = minDelay:frameDur:maxDelay; 

startTime = 0;
cnt = 1;

targetTimes = startTime + datasample(delays,1); 

while 1
    
    lastTime = targetTimes(end);
    
    newTargTime = lastTime + datasample(delays,1); % generate new target time
    
    if newTargTime > runDur-minDelay
        break % exit loop
    end
    
    % append new target time
    targetTimes = [targetTimes, newTargTime];
    
end

nTargs = length(targetTimes);