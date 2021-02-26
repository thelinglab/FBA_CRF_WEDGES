function [stair] = staircase_init(init_delta,step_size,ratio)
% initialize a staircase structure. 
%
% Inputs:
% init_delta: initial delta value
% step_size: size of step as a proportion (e.g. 0.02 is a 2% change)

stair.trialCnt = 1;
stair.acc = nan(1,2000);
stair.delta = nan(1,2000); % 2000 entries (shouldn't run out space)
stair.init_delta = init_delta; % initial value
stair.step_down = step_size; 
stair.step_up = step_size/ratio; % see Garcia-Perez 1998 Vis Research, pg 1871 for table
% the rule:
% delta_up = delta_down*(p/(1-p)) where p is the target accuracy
% but Garcia-Perez says this doesn't work that we

end

