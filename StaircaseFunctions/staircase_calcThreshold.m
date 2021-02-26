function [thresh_est,nReversals,nTrials] = staircase_calcThreshold(staircase_vals,reversals_criterion)
% inputs:
% staircase_vals = vector of staircase values
% reversals_criterion = number of reversals before calculating average
%
% outputs:
% thresh_est = threshold estimate
% reversal_cnt = how many reversals occured
% nTrials = number of trials after the nth reversal (where n is the
% reversals_criterion).

val = staircase_vals(~isnan(staircase_vals)); % remove nans

% preallocate
nEntries = length(val);
sign_change = nan(nEntries,1);
reversals = zeros(nEntries,1); 
aveIdx = zeros(nEntries,1);

% determine which trials the direction of the staircase changed
for e = 2:nEntries
    sign_change(e) = sign(val(e)-val(e-1));
end

% mark reversals
reversal_cnt = 0;
for e = 2:nEntries-1
   
    % mark trial to count in average if we've passed the reversal threshold
    if reversal_cnt >= reversals_criterion
        aveIdx(e) = 1;
    end
    
    if sign_change(e) == sign_change(e+1)
      reversals(e) = 0; 
   else
       reversals(e) = 1;
       reversal_cnt = reversal_cnt + 1; 
   end

end
aveIdx(end) = 1; % always include last value in average

% calculate outputs
nReversals = reversal_cnt;
nTrials = length(val(aveIdx == 1)); 
thresh_est = mean(val(aveIdx == 1)); % average across values after nth reversal

% for debugging
% plot(1:nEntries,val); hold on;
% scatter(1:nEntries,0.1*reversals,'r');
% scatter(1:nEntries,0.12*aveIdx,'b');