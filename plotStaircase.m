function plotStaircase(SubId)

reversal_thresh = 15;
col  = ['b','c','r','m'];

root = pwd;
StaircaseFunctionsDir = [root,'/StaircaseFunctions/'];
addpath(StaircaseFunctionsDir);

data_dir = ['Subject Data/',SubId,'/'];
fn = dir([data_dir,num2str(SubId),'_FBA_CRF_WEDGES_Staircase*.mat']); % grab the files for each subject
[nRuns,c]=size(fn);
fname = [data_dir,SubId,'_FBA_CRF_WEDGES_Staircase',num2str(nRuns),'.mat'];
load(fname,'staircase');

for cond = 1:4
    
    if cond == 1
        delta = staircase.upperVert.delta(~isnan(staircase.upperVert.delta));
    elseif cond == 2
        delta = staircase.lowerVert.delta(~isnan(staircase.lowerVert.delta));
    elseif cond == 3
        delta = staircase.upperHorz.delta(~isnan(staircase.upperHorz.delta));
    else
        delta= staircase.lowerHorz.delta(~isnan(staircase.lowerHorz.delta));
    end
    
    [thresh, reversals, trials_averaged] = staircase_calcThreshold(delta,reversal_thresh);
    
    % plot the staircase data
    y = delta;
    nTrials = length(y);
    t = 1:nTrials;
    plot(t,y,col(cond)); hold on;
    plot(nTrials-trials_averaged+1:nTrials,thresh*ones(1,trials_averaged),col(cond),'LineStyle','--');
    
end

xlabel('Trial')
ylabel('Change in SF');
legend('upper vert','','lower vert','','upper horz','','lower horz','');
