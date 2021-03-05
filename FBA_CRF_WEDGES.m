function FBA_CRF_WEDGES
% FBA_CRF (516 TRs, ~ 9 mins including before trigger is sent)
% Task: On each trial, the observer is cued to the upper or lower stimulus
% on the relevant side of the display. Report how many perturbations of
% spatial frequency occur in cued stimulus (0,1,2),
%
% Joshua J. Foster
% jjfoster@bu.edu
% Boston University


% To do:
% 1. add borders around the gratings on the relevant side so they're easily seggregated
% 2. The timing of the targets can be truly independent now
% given feedback at the end of trial? tricky because responses are quite
% slow
% 3. make textures half the visual field to speed up load times (lower
% priority)

 clear all

p.subID = 'S100';
p.runType = 0; % 0 = practice, 1 = actual run (just determines names of saved files)
p.runNum = 16;
p.surroundEccen = 2; % 1 or 2 % REVIEW
p.probeOri = 90; % 0 = vertical, 90 = horizontal
p.probeSide = 2; % 1 = probe on left, 2 = probe on right
p.targModulationRGB = 6; % modulation of gratings in RGB units
p.cueCB = 1; % 1 = pink = vertical, 2 = pink horizontal
p.eyeTrack = 0;
p.setup = 0; % 0 = my macbook, 1 = Display++, 2 = scanner at imaging center
echo off
KbName('UnifyKeyNames')
Screen('Preference', 'SkipSyncTests', 1);
input('hit enter to begin... ');

%% random number generator

% % seed the random generator
rng default % sets the generator to twister and sets the seed to 0 (as in matlab restarted)
rng shuffle % generates a new seed based on the clock
p.rngSettings = rng; % save to p-struct

%% setup file names

root = pwd;
data_dir = [root,'/Subject Data/',num2str(p.subID),'/'];

% subfolder for individual subject's data
if ~exist(data_dir,'dir')
    mkdir(data_dir);
end

% staircase run
if p.runType == 0
    fName = [data_dir,num2str(p.subID),'_FBA_CRF_WEDGES_Staircase',num2str(p.runNum),'.mat'];
end

% scanner run
if p.runType == 1
    fName = [data_dir,num2str(p.subID),'_FBA_CRF_WEDGES_Run',num2str(p.runNum),'.mat'];
end

% check that no such data file already exists
if exist(fName)
    Screen('CloseAll');
    msgbox('File already exists!', 'modal');
    return;
end

%% Setup display parameters
[keyboardIndices, productNames, ~] = GetKeyboardIndices;
if p.setup == 0 % My Macbook Pro
    w.whichScreen = 1;
    deviceString = 'Apple Internal Keyboard / Trackpad'; % macbook keyboard
    p.keyPressNumbers = [KbName('LeftArrow') KbName('RightArrow')];
    w.refreshRate = 60;
    w.screenWidth = 28.6; % my macbook
    w.vDist = 57; % viewing distance
    w.screenWidthPixels = 1440; % my macbook
    w.screenHeightPixels = 900;
    w.ScreenSizePixels = Screen('Rect', w.whichScreen);
end
if p.setup == 1 % Display++ for staircase session in Room 209
    w.whichScreen = 0;
    deviceString = 'Apple Keyboard'; % for both trigger and response box
    p.keyPressNumbers = [KbName('1') KbName('2')];
    w.refreshRate = 100;
    w.screenWidth = 52;      % screen width in cm
    w.vDist = 135;           % viewing distance
    w.screenWidthPixels = 1440;
    w.screenHeightPixels = 1080;
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    w.ScreenSizePixels = Screen('Rect', w.whichScreen);
end
if p.setup == 2 % Scanner
    w.whichScreen = 0;
    deviceString = '932'; % for both trigger and response box
    p.keyPressNumbers = [KbName('1!') KbName('2@')];
    w.refreshRate = 60;
    w.screenWidth = 42.7;   % screen width in cm (previously 43 cm)
    w.vDist = 99;           % viewing distance (previously 102 cm)
    w.screenWidthPixels = 1024;
    w.screenHeightPixels = 768; % Scanner display = [0 0 1024 768];
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    w.ScreenSizePixels = Screen('Rect', w.whichScreen);
end
for i=1:length(productNames)                                                % for each possible device
    if strcmp(productNames{i},deviceString)                                 % compare the name to the name you want
        deviceNumber=keyboardIndices(i);                                    % grab the correct id, and exit loop
        break;
    end
end
if deviceNumber==0 %error checking
    error('No device by that name was detected');
end

%% Screen parameters
p.pointSize = (2*atan2((2.54/72)/2, w.vDist))*(180/pi); % fontsize is in points (not pixel!) 1 letter point is 1/72 of an inch visual angle
w.screenVisAngle = 2*atan2d(w.screenWidth/2,w.vDist); % full with of monitor in degrees of visual angle
w.ppd = round(w.screenWidthPixels/w.screenVisAngle); % pixels per degree of visual angle
w.px = w.screenWidth/w.screenWidthPixels; % size of pixel in cm

% save current directory
root = pwd;
% add folder with called functions to path
CalledFunctionsDir = [root,'/CalledFunctions/'];
addpath(CalledFunctionsDir);
StaircaseFunctionsDir = [root,'/StaircaseFunctions/'];
addpath(StaircaseFunctionsDir);

%% General experiment parameters

% stimulus dimensions
p.fixSize = round(0.2*w.ppd);
p.SF = 1.5; % spatial frequency
p.gratingOuterAng = 70; % how far vertical boundaries are from horizontal meridian
p.gratingInnerAng = 0;  % how far horizontal boundaries are from horizontal meridiam
p.gratingInnerEdge = 1.0*w.ppd;
p.gratingOuterEdge = 8.5*w.ppd;
p.attendedInnerEdge = 1.0*w.ppd;
p.gapSize = 0.5*w.ppd; % size of gap between gratings on attended side
p.rolloffSD = 0.1*w.ppd;
p.eyeBoxSize = round(2.5*w.ppd);
p.fontSize = 40;
p.numPhaseSteps = 90;  % number of steps between 0 and 360 degrees, more increases load time

% orientations
p.oris = [0 90]; % 0 = vertical, 90 = horizontal
p.probeOri = 90; % probe is always horizontal

% Timing
p.refreshRate = w.refreshRate;
p.refreshCycle = 1/p.refreshRate;
p.frameRate = 10; % # of frames per second
p.frameDur = 1/p.frameRate;
p.adaptDur = 60; % duration of baseline adaptation
p.nAdaptFrames = p.adaptDur/p.frameDur;
p.nCueFrames = 5; % 5 = 500 ms
p.trialDur = 4;
p.nTrialFrames = p.trialDur/p.frameDur;
p.endBufferDur = 16;
p.eventSeqDur = 440;
p.runDur = p.adaptDur + p.eventSeqDur + p.endBufferDur;
p.nFrames = p.runDur/p.frameDur;
p.responseWindow = 2;
p.feedbackDur = 0.5; % for 500 ms at the end of the respose window

% timing of targets
p.targDur = 0.1;
p.targDur_frames = round(p.targDur/p.frameDur);
p.targLimits = [0.7 3.7]; % earliest and latest time grating can onset
p.targLimits_frames = round(p.targLimits./p.frameDur);
p.targ2Limits = [3.0 3.7]; % same for 2nd grating (must always occur toward end of trial)
p.targ2Limits_frames = round(p.targ2Limits./p.frameDur);
p.minSOA = 0.7; % min SOA between gratings
p.minSOA_frames = round(p.minSOA/p.frameDur);

% Color information
p.meanLum = 128;
p.targContrast = p.targModulationRGB/p.meanLum; 
p.bgCol = p.meanLum;
p.fixCol = 255;
p.txtCol = p.fixCol;
% cue colors both ~150 cd/m2 on projector in scanner
p.pink = [255,60,90];
p.blue = [40,175,240];
if  p.cueCB == 1
    p.upperCueCol = p.pink;
    p.lowerCueCol = p.blue;
end
if p.cueCB == 2
    p.upperCueCol = p.blue;
    p.lowerCueCol = p.pink;
end

% response keys
p.oneKey = p.keyPressNumbers(1);
p.twoKey = p.keyPressNumbers(2);

%% Generate contrast steps

stim.adaptContrast = 0.16;
stim.contrastSteps=[];
stim.octaves=[1 1.5 2 3];    % # of octaves of peaks/troughs from baseline

above = []; below = [];
for o = stim.octaves
    above = [above stim.adaptContrast*(2*o)];
    below = [stim.adaptContrast/(2*o) below];
end
stim.contrastSteps = [below stim.adaptContrast above];
stim.contrastLevels = [stim.contrastSteps];
p.nContrasts = length(stim.contrastLevels);
[m idx] = min(abs(stim.contrastLevels - stim.adaptContrast));
p.adaptContrastIdx = idx; 

%% Open window

AssertOpenGL;                                                               % make sure we're on an openGL-compatible machine
[window, p.sRect] = PsychImaging('OpenWindow', w.whichScreen, p.bgCol);
Screen('BlendFunction',window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);       % Turn on alpha blending. this makes drawing the stims much easier.
priorityLevel=MaxPriority(window);
Priority(priorityLevel);                                                    % set priority to max to prevent interruptions
HideCursor;

% force linear CLUT with ProPixx projector
w.OriginalCLUT = Screen('ReadNormalizedGammaTable', window);

if p.setup == 0
    load('MyMacbookPro_GammaCorrection.mat','calib'); % load gamma correction for my MBP
    w.linearizedCLUT = repmat(calib.gammaTable2,1,3); % save gamma correction to p struct
    Screen('LoadNormalizedGammaTable',window,w.linearizedCLUT);
else
    w.linearizedCLUT = linspace(0,1,256)'*ones(1,3);
    Screen('LoadNormalizedGammaTable',window,w.linearizedCLUT);
end

%% Create stimuli

% Create patches
p.xCenter = w.screenWidthPixels/2; p.yCenter = w.screenHeightPixels/2;
Screen('TextStyle', window, 1);
Screen('TextSize', window, p.fontSize);
bbox = Screen('TextBounds', window, 'Loading Stimuli... ');
newRect = CenterRectOnPoint(bbox, p.xCenter, p.yCenter);
Screen('DrawText', window, 'Loading Stimuli... ', newRect(1), newRect(2), p.txtCol);
Screen('Flip', window);

% Rects
p.bgRect = p.sRect; % background rect
p.eyeRect = CenterRectOnPoint([0 0 p.eyeBoxSize p.eyeBoxSize],p.xCenter,p.yCenter); % bounding box for eye position
p.fixRect = CenterRect([0 0 p.fixSize p.fixSize],p.bgRect);

%% make probe gratings

% set up phase values
lin=linspace(0,1,(p.numPhaseSteps+1));
p.phaseSteps=lin(1:(end-1))*2*pi;

p.stimSize = round2even(2*p.gratingOuterEdge + 4*p.rolloffSD); % allow 2x rolloffSD for each side
p.leftRect = CenterRectOnPoint([0 0 p.stimSize/2 p.stimSize+1], p.xCenter - p.stimSize/4, p.yCenter);
p.rightRect = CenterRectOnPoint([0 0 p.stimSize/2 p.stimSize+1], p.xCenter + p.stimSize/4, p.yCenter);
if p.probeSide == 1
    p.probeRect = p.leftRect;
    p.attRect = p.rightRect;
else
    p.probeRect = p.rightRect;
    p.attRect = p.leftRect;
end

[Xc,Yc] = meshgrid(1:(p.stimSize/2), (-p.stimSize/2):(p.stimSize/2));
eccen = sqrt((Yc).^2+(Xc).^2);
ang = nan(size(Xc));
for y = 1:p.stimSize+1
    for x = 1:p.stimSize/2
        ang(y,x) = rad2deg(cart2pol(Xc(y,x),Yc(y,x)));
    end
end

% stim on unattended side
probe_aperture = zeros(size(Xc));
probe_aperture(eccen >= p.gratingInnerEdge & eccen <= p.gratingOuterEdge & abs(ang) <= p.gratingOuterAng) = 1; % hemifield
if p.probeSide == 1
    probe_aperture = fliplr(probe_aperture);
end
probe_aperture = imgaussfilt(probe_aperture,p.rolloffSD);

% upper aperture
upper_aperture = zeros(size(Xc));
upper_aperture(eccen >= p.gratingInnerEdge & eccen <= p.gratingOuterEdge & Yc < -p.gapSize/2 & ang > -p.gratingOuterAng) = 1;
if p.probeSide == 2
    upper_aperture = fliplr(upper_aperture);
end
upper_aperture = imgaussfilt(upper_aperture,p.rolloffSD);

% lower aperture
lower_aperture = zeros(size(Xc));
lower_aperture(eccen >= p.gratingInnerEdge & eccen <= p.gratingOuterEdge & Yc > p.gapSize/2 & ang < p.gratingOuterAng) = 1;
if p.probeSide == 2
    lower_aperture = fliplr(lower_aperture);
end
lower_aperture = imgaussfilt(lower_aperture,p.rolloffSD);

tic
for s=1:p.numPhaseSteps
    
    sinusoid = cos(2*pi*p.SF/w.ppd*(Xc.*cos(p.probeOri*(pi/180))+Yc.*sin(p.probeOri*(pi/180)))-p.phaseSteps(s));
    
    for c = 1:p.nContrasts
        
        tmp = nan(size(Xc,1),size(Xc,2),2);
        tmp(:,:,1) = sinusoid.*p.meanLum(1)*stim.contrastLevels(c)+p.meanLum(1);
        tmp(:,:,2) = probe_aperture*255;
        probeText{c,s} = Screen('MakeTexture',window,tmp);
        
    end
    
    for ori = 1:2
        
        sinusoid = cos(2*pi*p.SF/w.ppd*(Xc.*cos(p.oris(ori)*(pi/180))+Yc.*sin(p.oris(ori)*(pi/180)))-p.phaseSteps(s));
        
        tmp = nan(size(Xc,1),size(Xc,2),2);
        tmp(:,:,1) = sinusoid.*p.meanLum(1)*stim.contrastLevels(5)+p.meanLum(1);
        tmp(:,:,2) = upper_aperture*255;
        upperText{ori,s} = Screen('MakeTexture',window,tmp);
        
        tmp = nan(size(Xc,1),size(Xc,2),2);
        tmp(:,:,1) = sinusoid.*p.meanLum(1)*stim.contrastLevels(5)+p.meanLum(1);
        tmp(:,:,2) = lower_aperture*255;
        lowerText{ori,s} = Screen('MakeTexture',window,tmp);
        
    end
    
end

toc


%% trial sequence

stim.adaptFlipTimes = 0:p.frameDur:p.adaptDur-p.frameDur;
stim.nAdaptFrames = length(stim.adaptFlipTimes);

% load paradigm file
par_dir = [data_dir,'ParFiles/'];
runNumString = num2str(p.runNum); l = length(runNumString);
if l == 1
    parFileName = ['PAR-00',runNumString,'.par'];
elseif l == 2
    parFileName = ['PAR-0',runNumString,'.par'];
else
    error('cannot find par file by that name');
end
p.paradigmFile = [par_dir,parFileName];
fid = fopen(p.paradigmFile);
p.paradigmStruct = textscan(fid,'%f %f %f %f %s', 'Delimiter', '\t'); fclose('all');
% check that length of run in par file matches p.eventSeqDur
eventSeqDur = p.paradigmStruct{1}(end)+p.paradigmStruct{3}(end);
if p.eventSeqDur ~= eventSeqDur
    error('p.eventSeqDur does not match par file')
end
stim.absTimes = p.paradigmStruct{1}; % cumulative time
stim.eventDur = p.paradigmStruct{3}; % duration of each event
stim.condNumber = p.paradigmStruct{2}; % conditions number
stim.condName = p.paradigmStruct{5}; % condition

% duration of null/topup period
stim.topupDur = stim.eventDur(stim.condNumber == 0);

% remove the nulls and blank trials
stim.absTimes = stim.absTimes(stim.condNumber ~= 0);
stim.condName = stim.condName(stim.condNumber ~= 0);
stim.condNumber = stim.condNumber(stim.condNumber ~= 0);
stim.nTrials = length(stim.absTimes);

% calculate number of frames in each null/topup period
stim.nTopupFrames = nan(stim.nTrials,1);
for t = 1:stim.nTrials
    stim.nTopupFrames(t) = stim.topupDur(t)/p.frameDur;
end

% determine direction of each probe
stim.cuedOri = nan(stim.nTrials,1);
stim.uncuedOri = nan(stim.nTrials,1);
% vertical probe
if p.probeOri == 0
    stim.cuedOri(stim.condNumber <= 10) = 0; % match (attended  = vertical)
    stim.uncuedOri(stim.condNumber <= 10) = 90; 
    stim.cuedOri(stim.condNumber > 10) = 90;  % mismatch (attended = horizontal)
    stim.uncuedOri(stim.condNumber > 10) = 0; 
end
% horizontal probe
if p.probeOri == 90
    stim.cuedOri(stim.condNumber <= 10) = 90; % match (attended = horizontal)
    stim.uncuedOri(stim.condNumber <= 10) = 0;
    stim.cuedOri(stim.condNumber > 10) = 0;   % mismatch (attended  = vertical)
    stim.uncuedOri(stim.condNumber > 10) = 90; 
end

% determine position of attended stimulus (upper or lower)
stim.cuedPos = nan(stim.nTrials,1);
for cond = 1:max(stim.condNumber)
    condIdx = find(stim.condNumber == cond);
    stim.cuedPos(condIdx) = Shuffle([1 2]); % a trial with each position attended for each condition
end

% determine orientation of each grating on atended side
stim.upperOri = nan(stim.nTrials,1);
stim.upperOri(stim.cuedPos == 1) = stim.cuedOri(stim.cuedPos == 1);
stim.upperOri(stim.cuedPos == 2) = stim.uncuedOri(stim.cuedPos == 2);

stim.lowerOri = nan(stim.nTrials,1);
stim.lowerOri(stim.cuedPos == 2) = stim.cuedOri(stim.cuedPos == 2);
stim.lowerOri(stim.cuedPos == 1) = stim.uncuedOri(stim.cuedPos == 1);

% determine probe contrast

stim.probeContrastIdx = nan(stim.nTrials,1);
for t = 1:stim.nTrials
    stim.probeContrastIdx(t) = str2num(stim.condName{t}(end)); % grabbing contrast level from stim.condName struct
end

% determine trial times
stim.trialTimes = p.adaptDur + stim.absTimes;

% 0, 1, or 2 targs equally often for each orientation, randomized independently 

% number of vert targets on each trial
tmp = [zeros(1,13),ones(1,13),2*ones(1,13),datasample(0:2,1)]; % 13 of each, plus one extra randomly determined
permIdx = randperm(stim.nTrials);
stim.nVertTargs = tmp(permIdx);
% number of horz targets on each trial
tmp = [zeros(1,13),ones(1,13),2*ones(1,13),datasample(0:2,1)]; % 13 of each, plus one extra randomly determined
permIdx = randperm(stim.nTrials);
stim.nHorzTargs = tmp(permIdx);

stim.nTargs = nan(1,stim.nTrials); 
stim.nTargs(stim.cuedOri == 0) = stim.nVertTargs(stim.cuedOri == 0);
stim.nTargs(stim.cuedOri == 90) = stim.nHorzTargs(stim.cuedOri == 90); 

% preallocate vectors to save response information
stim.resp = zeros(1,stim.nTrials);
stim.acc = zeros(1,stim.nTrials);
stim.rt = nan(1,stim.nTrials);

% determine when gratings are presented
stim.vertTarg1On = zeros(stim.nTrials,p.nTrialFrames);
stim.vertTarg2On = zeros(stim.nTrials,p.nTrialFrames);
stim.horzTarg1On = zeros(stim.nTrials,p.nTrialFrames);
stim.horzTarg2On = zeros(stim.nTrials,p.nTrialFrames);
stim.vertTarg1Frame = zeros(stim.nTrials,p.nTrialFrames);
stim.vertTarg2Frame = zeros(stim.nTrials,p.nTrialFrames);
stim.horzTarg1Frame = zeros(stim.nTrials,p.nTrialFrames);
stim.horzTarg2Frame = zeros(stim.nTrials,p.nTrialFrames);

for t = 1:stim.nTrials
        
    % randomly sample phase for each grating
    stim.vertGratingPhase{t} = datasample(p.phaseSteps,2); % always generating 2 phase values, but only using both on trials with two gratings
    stim.horzGratingPhase{t} = datasample(p.phaseSteps,2);
    
    
    vertOnsets = []; horzOnsets = []; 
    while 1

        % no vert target
        if stim.nVertTargs(t) == 0
        vertOnsets = [];
        end
        
        % one vert target
        if stim.nVertTargs(t) == 1;
        vertOnsets = datasample(p.targLimits_frames(1):p.targLimits_frames(2),1);
        end
        
        % two vert targets
        if stim.nVertTargs(t) == 2
        targ2onset = datasample(p.targ2Limits_frames(1):p.targ2Limits_frames(2),1);
        targ1onset = datasample(p.targLimits_frames(1):targ2onset-p.minSOA_frames,1);
        vertOnsets = [targ1onset targ2onset];
        end
        
        % no horz target
        if stim.nHorzTargs(t) == 0
            horzOnsets = [];
        end
        
        % one horz target
        if stim.nHorzTargs(t) == 1
            horzOnsets = datasample(p.targLimits_frames(1):p.targLimits_frames(2),1);
        end
        
        % two horz targets
        if stim.nHorzTargs(t) == 2
            targ2onset = datasample(p.targ2Limits_frames(1):p.targ2Limits_frames(2),1);
            targ1onset = datasample(p.targLimits_frames(1):targ2onset-p.minSOA_frames,1);
            horzOnsets = [targ1onset targ2onset];
        end
        
        onsets = sort([vertOnsets, horzOnsets]);
        
        if length(onsets) < 2
            minSOA = p.minSOA+100; % make minSOA bigger than p.minSOA if there is zero or one onsets
        else
            minSOA = calcMinSOA(onsets);
        end
        
        if minSOA > p.minSOA_frames
            break
        end
        
    end
    
    stim.vertTargOnsets{t} = vertOnsets;
    stim.horzTargOnsets{t} = horzOnsets;
    stim.onsets{t} = onsets;
    
    
    if stim.nVertTargs(t) == 1
        stim.vertTarg1On(t,stim.vertTargOnsets{t}(1):stim.vertTargOnsets{t}(1)+p.targDur_frames-1) = 1;
        stim.vertTarg1Frame(t,stim.vertTargOnsets{t}(1):stim.vertTargOnsets{t}(1)+p.targDur_frames-1) = 1:p.targDur_frames;
    end
    
    if stim.nVertTargs(t) == 2
        stim.vertTarg1On(t,stim.vertTargOnsets{t}(1):stim.vertTargOnsets{t}(1)+p.targDur_frames-1) = 1;
        stim.vertTarg1Frame(t,stim.vertTargOnsets{t}(1):stim.vertTargOnsets{t}(1)+p.targDur_frames-1) = 1:p.targDur_frames;
        stim.vertTarg2On(t,stim.vertTargOnsets{t}(2):stim.vertTargOnsets{t}(2)+p.targDur_frames-1) = 1;
        stim.vertTarg2Frame(t,stim.vertTargOnsets{t}(2):stim.vertTargOnsets{t}(2)+p.targDur_frames-1) = 1:p.targDur_frames;
    end
    
    if stim.nHorzTargs(t) == 1
        stim.horzTarg1On(t,stim.horzTargOnsets{t}(1):stim.horzTargOnsets{t}(1)+p.targDur_frames-1) = 1;
        stim.horzTarg1Frame(t,stim.horzTargOnsets{t}(1):stim.horzTargOnsets{t}(1)+p.targDur_frames-1) = 1:p.targDur_frames;
    end
    
    if stim.nHorzTargs(t) == 2
        stim.horzTarg1On(t,stim.horzTargOnsets{t}(1):stim.horzTargOnsets{t}(1)+p.targDur_frames-1) = 1;
        stim.horzTarg1Frame(t,stim.horzTargOnsets{t}(1):stim.horzTargOnsets{t}(1)+p.targDur_frames-1) = 1:p.targDur_frames;
        stim.horzTarg2On(t,stim.horzTargOnsets{t}(2):stim.horzTargOnsets{t}(2)+p.targDur_frames-1) = 1;
        stim.horzTarg2Frame(t,stim.horzTargOnsets{t}(2):stim.horzTargOnsets{t}(2)+p.targDur_frames-1) = 1:p.targDur_frames;
    end
    
end
        
%% Eye link setup

if p.eyeTrack == 1
    
    if p.runType == 0
        edfFileFinal = [num2str(p.subID),'_FBA_CRF_WEDGES_Staircase',num2str(p.runNum),'.edf'];
    end
    if p.runType == 1
        edfFileFinal = [num2str(p.subID),'_FBA_CRF_WEDGES_Run',num2str(p.runNum),'.edf'];
    end
    
    edf_filename = 'tmp.edf';
    
    [el] = eyeTrackingOn_500Hz(window, edf_filename, p.bgRect, p.eyeRect, w.ppd, p.bgCol, p.txtCol);
    Screen('Flip', window);
end

%% wait for trigger
Screen('FillRect',window,p.bgCol,p.bgRect);
Screen('FillOval',window,p.fixCol,p.fixRect);
DrawFormattedText(window, '~', w.ppd, w.ppd, p.txtCol, p.bgCol);
Screen('Flip', window);

% Create a queue which records keypresses from button box only
triggerKey = KbName('5%');
keylist = zeros(1,256); % keys for KbQueueCreate
keylist(p.keyPressNumbers) = 1; % only monitor p.keyPressNumbers

fprintf('Waiting for trigger... \n')

% Wait for scanner trigger or mouse click
if p.setup == 2 % only wait for trigger if we're running in the scanner
    KbTriggerWait(triggerKey, deviceNumber);
else
    GetClicks;
end
fprintf('Trigger detected! \n')
PsychHID('KbQueueCreate', deviceNumber, keylist);
PsychHID('KbQueueStart', deviceNumber);

% start eye tracker recording
if p.eyeTrack
    Eyelink('StartRecording');
    % mark zero-plot time in EDF file
    Eyelink('Message', 'RunStart');  % may be a lag in recording (not a problem because the baseline screen is much longer than needed)
end

% determine event timing as soon as the trigger pulse is received
presTime.runStart = GetSecs; time.runStart = presTime.runStart;
presTime.runEnd = presTime.runStart + p.runDur;

presTime.adaptFlipTimes = presTime.runStart + stim.adaptFlipTimes;
presTime.trialStart = presTime.runStart + stim.trialTimes;
presTime.topupStart = presTime.trialStart + p.trialDur;
for t = 1:stim.nTrials
    presTime.trialFlipTimes{t} = presTime.trialStart(t):p.frameDur:presTime.trialStart(t)+p.trialDur-p.frameDur;
    presTime.topupFlipTimes{t} = presTime.topupStart(t):p.frameDur:presTime.topupStart(t)+stim.topupDur(t)-p.frameDur;
end
endBufferStart = presTime.runStart + p.adaptDur + p.eventSeqDur;
presTime.endBufferFlipTimes = endBufferStart:p.frameDur:presTime.runEnd-p.frameDur;
stim.nEndBufferFrames = length(presTime.endBufferFlipTimes);

% variables for recording psychtoolbox timing
time.adapt.vblstamp = nan(stim.nAdaptFrames,1);
time.adapt.onset =  nan(stim.nAdaptFrames,1);
time.adapt.flipstamp = nan(stim.nAdaptFrames,1);
time.adapt.missed = nan(stim.nAdaptFrames,1);

time.trial.vblstamp = nan(stim.nTrials,p.nTrialFrames);
time.trial.onset =  nan(stim.nTrials,p.nTrialFrames);
time.trial.flipstamp = nan(stim.nTrials,p.nTrialFrames);
time.trial.missed = nan(stim.nTrials,p.nTrialFrames);

time.topupAdapt.vblstamp = {};
time.topupAdapt.onset =  {};
time.topupAdapt.flipstamp = {};
time.topupAdapt.missed = {};

time.endBuffer.vblstamp = nan(stim.nEndBufferFrames,1);
time.endBuffer.onset =  nan(stim.nEndBufferFrames,1);
time.endBuffer.flipstamp = nan(stim.nEndBufferFrames,1);
time.endBuffer.missed = nan(stim.nEndBufferFrames,1);

%% begin task

%% base adaptation period
cueFrames = p.nAdaptFrames-p.nCueFrames:p.nAdaptFrames;

for f = 1:p.nAdaptFrames
    
    %% draw display
    
    % background
    Screen('FillRect',window,p.bgCol,p.bgRect);
    
    % determine phase of gratings
     currPhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
        
    % draw surround grating
    Screen('DrawTexture',window, probeText{p.adaptContrastIdx,currPhase(1)},[],p.stimRect,0);
    
    % fixation
    if ismember(f,cueFrames)
        if stim.cuedPos(1) == 1
            cueCol = p.upperCueCol;
        end
        if stim.cuedPos(1) == 2
            cueCol = p.lowerCueCol;
        end
        Screen('FillOval',window,cueCol,p.fixRect);
    else
        Screen('FillOval',window,p.fixCol,p.fixRect);
    end
    
    Screen('DrawingFinished',window);
    
    [time.adapt.vblstamp(f), time.adapt.onset(f), time.adapt.flipstamp(f), time.adapt.missed(f)] = Screen('Flip',window,presTime.adaptFlipTimes(f));
    
    
end

%% trial loop

for t = 1:stim.nTrials
    
    if stim.probeContrastIdx(t) == 0
        contrast = 0;
    else
        contrast = 100*stim.contrastLevels(stim.probeContrastIdx(t));
    end
    fprintf('Trial = %.0f ,',t);
    fprintf(stim.condName{t});
    fprintf(' , Cued position = %.0f, Targets = %.0f, Contrast = %.1f, ',stim.cuedPos(t),stim.nTargs(t),contrast)
    trialFlipTimes = presTime.trialFlipTimes{t};
    topupFlipTimes = presTime.topupFlipTimes{t};
    PsychHID('KbQueueFlush',deviceNumber); % flush any response made before trial began
    
    %% make target textures for trial
    
    upperTarg{2} = [];
    lowerTarg{2} = [];
    probeTarg{2} = [];
    
    tic
    
    for targNum = 1:2
        
        randPhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
        sinusoid = cos(2*pi*(p.SF - stim.upperDelta(t))/w.ppd*(Xc.*cos(stim.upperOri(t)*(pi/180))+Yc.*sin(stim.upperOri(t)*(pi/180)))-randPhase);
        tmp = nan(size(Xc,1),size(Xc,2),2);
        tmp(:,:,1) = sinusoid.*p.meanLum(1)*stim.contrastLevels(5)+p.meanLum(1);
        tmp(:,:,2) = upper_aperture*255;
        upperTarg{targNum} = Screen('MakeTexture',window,tmp);
        stim.upperTargPhase(t,targNum) = randPhase;
        
        randPhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
        sinusoid = cos(2*pi*(p.SF - stim.lowerDelta(t))/w.ppd*(Xc.*cos(stim.lowerOri(t)*(pi/180))+Yc.*sin(stim.lowerOri(t)*(pi/180)))-randPhase);
        tmp = nan(size(Xc,1),size(Xc,2),2);
        tmp(:,:,1) = sinusoid.*p.meanLum(1)*stim.contrastLevels(5)+p.meanLum(1);
        tmp(:,:,2) = lower_aperture*255;
        lowerTarg{targNum} = Screen('MakeTexture',window,tmp);
        stim.lowerTargPhase(t,targNum) = randPhase;
        
        randPhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
        sinusoid = cos(2*pi*(p.SF + stim.probeDelta(t))/w.ppd*(Xc.*cos(p.probeOri*(pi/180))+Yc.*sin(p.probeOri*(pi/180)))-randPhase);
        tmp = nan(size(Xc,1),size(Xc,2),2);
        tmp(:,:,1) = sinusoid.*p.meanLum(1)*contrast+p.meanLum(1);
        tmp(:,:,2) = probe_aperture*255;
        probeTarg{targNum} = Screen('MakeTexture',window,tmp);
        stim.probePhase(t,targNum) = randPhase;
        
    end
    
    time.makeTextures(s) = toc; % save time, so I can check that it didn't disrupt flip times
    
    
    %% trial period
    
    for f = 1:p.nTrialFrames
        
        % background
        Screen('FillRect',window,p.bgCol,p.bgRect);
        
        % determine phase of gratings
        currPhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
        unattendedPhase = round((rand(1,1)*(p.numPhaseSteps-1))+1); % REVIEW
                       
         % draw surround grating
         if stim.probeContrastIdx(t) ~= 0
             Screen('DrawTexture',window, probeText{stim.probeContrastIdx(t),currPhase(1)},[],p.stimRect,0);
         end
         if stim.upperOri(t) == 0
             upperOri = 1;
             lowerOri = 2;
             upperPhase = unattendedPhase;
             lowerPhase = currPhase;
         end
         if stim.upperOri(t) == 90
             upperOri = 2;
             lowerOri = 1;
             upperPhase = currPhase;
             lowerPhase = unattendedPhase;
         end
         Screen('DrawTexture',window, lowerText{lowerOri,lowerPhase},[],p.stimRect,0);
         Screen('DrawTexture',window, upperText{upperOri,upperPhase},[],p.stimRect,0);
         
        % fixation dot
        if stim.cuedPos(t) == 1
            cueCol = p.upperCueCol;
        end
        if stim.cuedPos(t) == 2
            cueCol = p.lowerCueCol;
        end
        Screen('FillOval',window,cueCol,p.fixRect);
        
        Screen('DrawingFinished',window);
        
        % flip display
        [time.trial.vblstamp(t,f), time.trial.onset(t,f), time.trial.flipstamp(t,f), time.trial.missed(t,f)] = Screen('Flip',window,trialFlipTimes(f));
        
        if f == 1
            % send event marker
            if p.eyeTrack
                Eyelink('message','TrialStart');
            end
        end
        
    end
    
    %% topup adaptation
    
    respStartTime = GetSecs;
    responseDeadline = respStartTime + p.responseWindow;
    cueFrames = stim.nTopupFrames(t)-p.nCueFrames:stim.nTopupFrames(t);
    
    for f = 1:stim.nTopupFrames(t)
        
        % background
        Screen('FillRect',window,p.bgCol,p.bgRect);
        
        % probe stimulus
            currPhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
      
        % draw surround grating
        Screen('DrawTexture',window, probeText{p.adaptContrastIdx,currPhase(1)},[],p.stimRect,0);
        
        % fixation dot
        if ismember(f,cueFrames) && t < stim.nTrials
            if stim.cuedPos(t+1) == 1
                cueCol = p.upperCueCol;
            end
            if stim.cuedPos(t+1) == 2
                cueCol = p.lowerCueCol;
            end
            Screen('FillOval',window,cueCol,p.fixRect);
        else
            Screen('FillOval',window,p.fixCol,p.fixRect);
        end
        
        Screen('DrawingFinished',window);
        
        % flip display
        [time.topupAdapt.vblstamp{t}(f), time.topupAdapt.onset{t}(f), time.topupAdapt.flipstamp{t}(f), time.topupAdapt.missed{t}(f)] = Screen('Flip',window,topupFlipTimes(f));
        
        if f == 1
            % send event marker
            if p.eyeTrack
                Eyelink('message','TopupStart');
            end
        end
        
        %% monitor for response
        
        while 1
            
            % break loop after frame duration
            if GetSecs >= topupFlipTimes(f)+p.frameDur
                break
            end
            
            % if responseDeadline exceeded....
            if GetSecs > responseDeadline
                % determine accuracy
                if stim.resp(t) == stim.nTargs(t)
                    stim.acc(t) = 1;
                end
                fprintf(' Resp = %0.f, Acc = %.0f\n',stim.resp(t),stim.acc(t));
                responseDeadline = nan; % set responseDeadline to a nan
            end
            
            % check for response if haven't hit the response deadline
            if GetSecs <= responseDeadline
                
                [pressed, keyIndex] = PsychHID('KbQueueCheck', deviceNumber);
                
                if pressed
                    
                    rt = GetSecs - respStartTime; % calculate rt
                    
                    [m, keyCode] = max(keyIndex);
                    
                    PsychHID('KbQueueFlush',deviceNumber);
                    
                    % save response time
                    stim.rt(t) = rt;
                    
                    % save response
                    if keyCode == p.oneKey
                        stim.resp(t) = 1;
                    end
                    if keyCode == p.twoKey
                        stim.resp(t) = 2;
                    end
                    
                end
                
            end
            
        end

    end
     
end

%% final null period

for f = 1:stim.nEndBufferFrames
    
    % background
    Screen('FillRect',window,p.bgCol,p.bgRect);
    
    % determine phase of gratings
    currPhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
    
    % draw surround grating
    Screen('DrawTexture',window, probeText{p.adaptContrastIdx,currPhase(1,1)},[],p.stimRect,0);
    
    % fixation
    Screen('FillOval',window,p.fixCol,p.fixRect);
    
    Screen('DrawingFinished',window);
    
    % flip display
    [time.endNull.vblstamp(f), time.endNull.onset(f), time.endNull.flipstamp(f), time.endNull.missed(f)] = Screen('Flip',window,presTime.endBufferFlipTimes(f));
    
end

WaitSecs('UntilTime',presTime.runEnd);
time.runEnd = GetSecs; % save end of run time

%% Wakefulness rating

PsychHID('KbQueueFlush', deviceNumber);
Screen('TextSize', window, 20);
xCenter = w.ScreenSizePixels(3)/2; yCenter = w.ScreenSizePixels(4)/2;
startloc = 2;
targLocations = linspace(-xCenter/2,xCenter/2,4);
targText = [{'Drowsy/asleep'} {'Somewhat drowsy'} {'Mostly awake'} {'Fully awake'}];
Text = 'How awake were you during this scan?';
startRatingTime = GetSecs;
while 1
    if GetSecs > startRatingTime+4
        break;
    end
    [pressed, firstpress] = PsychHID('KbQueueCheck', deviceNumber);
    PsychHID('KbQueueFlush',deviceNumber);
    [~, indice]=max(firstpress);
    if pressed
        if p.keyPressNumbers(1) == indice
            startloc = startloc - 1;
            if startloc < 1
                startloc = 1;
            end
        else
            startloc = startloc + 1;
            if startloc > 4
                startloc = 4;
            end
        end
    end
    DrawFormattedText(window, Text, 'center', yCenter-w.ppd*2, p.txtCol);
    for n = 1:numel(targLocations)
        textRect = CenterRectOnPoint(Screen('TextBounds', window, targText{n}), xCenter-targLocations(n), yCenter);
        DrawFormattedText(window, targText{n}, textRect(1), yCenter-w.ppd*1, p.txtCol);
        Screen('FillOval', window, p.txtCol, CenterRectOnPoint([0 0 round(w.ppd) round(w.ppd)], xCenter+targLocations(n), yCenter+w.ppd/2));
        Screen('FillOval', window, p.bgCol, CenterRectOnPoint([0 0 round(w.ppd)-4 round(w.ppd)-4], xCenter+targLocations(n), yCenter+w.ppd/2));
        numRect = CenterRectOnPoint(Screen('TextBounds', window, num2str(n)), xCenter+targLocations(n), yCenter+w.ppd/2);
        DrawFormattedText(window, num2str(n), numRect(1), numRect(2), p.txtCol);
    end
    Screen('FillOval', window, [255 0 0], CenterRectOnPoint([0 0 round(w.ppd) round(w.ppd)], xCenter+targLocations(startloc), yCenter+w.ppd/2));
    Screen('FillOval', window, p.bgCol, CenterRectOnPoint([0 0 round(w.ppd)-4 round(w.ppd)-4], xCenter+targLocations(startloc), yCenter+w.ppd/2));
    numRect = CenterRectOnPoint(Screen('TextBounds', window, num2str(startloc)), xCenter+targLocations(startloc), yCenter+w.ppd/2);
    DrawFormattedText(window, num2str(startloc), numRect(1), numRect(2), p.txtCol);
    Screen('Flip', window);
end

wake.WakefullnessRating = startloc;
wake.WakefullnessLabels = fliplr(targText);

%% Clean up

% save the data
save(fName,'p','w','stim','presTime','time','staircase','wake');

if p.eyeTrack == 1
    % reset so tracker uses default calibration for other experiments
    Eyelink('command', 'generate_default_targets = yes')
    Eyelink('Command', 'set_idle_mode');
    Eyelink('command', ['screen_pixel_coords' num2str(w.ScreenSizePixels)])
    Eyelink('StopRecording');
    Eyelink('CloseFile');
    transferEDF(edf_filename,edfFileFinal)
end

% restore clut
Screen('LoadNormalizedGammaTable', window, w.OriginalCLUT);
Screen('CloseAll');

PsychHID('KbQueueStop', deviceNumber);
ShowCursor;
Priority(0);

% print delta for the run
fprintf('Accuracy =  %.2f \n',100*mean(stim.acc));