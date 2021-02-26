function [stair] = staircase_update(stair,delta_sugg,acc)
% add accuracy for current trial to the structure, and advance trial
% counter
% inputs:
% staircase structure
% delta_sugg is suggested delta (not rounded)
% accuracy = 0 is incorrect, 1 is correct

t = stair.trialCnt; % get current trial count
stair.acc(t) = acc; % add accuracy
stair.delta(t) = delta_sugg;
stair.trialCnt = t+1;  % advance trial count


end






