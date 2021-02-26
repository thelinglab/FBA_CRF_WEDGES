function [delta] = staircase_getDelta(stair,maxChange)
% access staircase structure to work out the next value to present. 
% Inputs:
% stair: the staircase structure
% maxChange: biggest possible change
% Ouputs:
% [delta for trial]

t = stair.trialCnt; % get trial count

if t == 1 % for first trial use initial value
    delta = stair.init_delta;
end

if t > 1  % otherwise, apply the staircase algorithm
    
    lastDelta = stair.delta(t-1);
    lastAcc = stair.acc(t-1);

    % if last trial incorrect, make change bigger
    if lastAcc == 0
        delta = lastDelta*(1+stair.step_up); % go up by 
    end
    
    % if last trial correct, make change smaller
    if lastAcc == 1
        delta = lastDelta*(1-stair.step_down); 
    end     
end

% present biggest possible change if staircase is "out of bounds"
if delta > maxChange
    delta = maxChange;
end