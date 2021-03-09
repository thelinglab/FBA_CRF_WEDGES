function FBA_CRF_WEDGES_Localizer
% Localizer scan (272 TRs, ~5 mins including pre-scan time)
% Task: press any key when fixation dot dims.
%
% Joshua J. Foster
% jjfoster@bu.edu
% Boston University

% general settings
clear
p.subID = 'S100';
p.scan = 0; % 0 = running on laptop, 1 = running in scanner
p.runType = 0; % 0 = demo, 1 = actual data acquisition (determines name of saved file)
p.runNum = 1; 
p.eyeTrack = 0;
p.probeSide = 1; % 1 = irrelevant probe on left, 2 = probe on right

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

% demo run
if p.runType == 0
    fName = [data_dir,num2str(p.subID),'_FBA_CRF_WEDGES_Localizer_Demo',num2str(p.runNum),'.mat'];
end

% scanner run
if p.runType == 1
    fName = [data_dir,num2str(p.subID),'_FBA_CRF_WEDGES_Localizer_Run',num2str(p.runNum),'.mat'];
end

% check that no such data file already exists
if exist(fName)
    Screen('CloseAll');
    msgbox('File already exists!', 'modal');
    return;
end

%% Setup display parameters
[keyboardIndices, productNames, ~] = GetKeyboardIndices;
if p.scan == 0
    w.whichScreen = 1; 
    deviceString = 'Apple Internal Keyboard / Trackpad'; % macbook keyboard
    p.keyPressNumbers = [KbName('1!') KbName('2@')]; 
    w.refreshRate = 60;
    w.screenWidth = 28.6; % scanner macbook
    w.vDist = 57; % viewing distance
    w.screenWidthPixels = 1400;
    w.screenHeightPixels = 900; % scanner macbook
    w.ScreenSizePixels = Screen('Rect', w.whichScreen);                        
else
    w.whichScreen = 0;
    deviceString = '932'; % for both trigger and response box 
    p.keyPressNumbers =[KbName('1!') KbName('2@')]; 
    w.refreshRate = 60;
    w.screenWidth = 42.7;   % screen width in cm (previously 43 cm)
    w.vDist = 99;           % viewing distance (previously 102 cm)
    w.screenWidthPixels = 1024;
    w.screenHeightPixels = 768; % Scanner display = [0 0 1024 768];
%     Datapixx('Open'); 
%     Datapixx('RegWrRd');
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible')
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
p.gratingInnerEdge = 1.0*w.ppd;
p.gratingOuterEdge = 8.5*w.ppd;
p.gapSize = 0.5*w.ppd; % size of gap between gratings on attended side
p.rolloffSD = 0.1*w.ppd;
p.eyeBoxSize = round(2.5*w.ppd);
p.fontSize = 40;
p.numPhaseSteps = 90; % number of steps between 0 and 360 degrees, more increases load time
p.ori = 90;


% Timing
p.refreshRate = 60;
p.refreshCycle = 1/p.refreshRate;
p.stimFreq = 10; % stimulation frequency
p.frameDur = 1/p.stimFreq;
p.blockDur = 16; 
p.blocksPerStim = 8;
p.cond = 2;
p.nBlocks = p.cond*p.blocksPerStim+1; % +1 for final baseline
p.runDur = p.blockDur*p.nBlocks; 
p.nFrames = p.runDur/p.frameDur; % total stimulus frames in run

% timing of fixation targets
p.minTargDelay = 3; % minimum target-to-target SOA 
p.maxTargDelay = 5; % maximum target-to-target SOA
p.responseWindow = 2; % time to respond to target

% Color information
p.meanLum = 128;
p.bgCol = p.meanLum;
p.maxContrastRGB = 127; % maximum modulation in RGB units (corresponds to 1-255 RGB units)
p.white = 255;
p.txtCol = p.white;
p.fixTargCol = 190; % color of fixation targets

%% Open window

AssertOpenGL;                                                               % make sure we're on an openGL-compatible machine
[window, p.sRect] = PsychImaging('OpenWindow', w.whichScreen, p.bgCol);
Screen('BlendFunction',window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);       % Turn on alpha blending. this makes drawing the stims much easier.
priorityLevel=MaxPriority(window);
Priority(priorityLevel);                                                    % set priority to max to prevent interruptions
HideCursor;

% force linear CLUT with ProPixx projector
w.OriginalCLUT = Screen('ReadNormalizedGammaTable', window);
if p.scan == 1
    Screen('LoadNormalizedGammaTable', window, linspace(0,1,256)'*ones(1,3));
end

%% Create stimuli

% Create patches
xCenter = w.screenWidthPixels/2; yCenter = w.screenHeightPixels/2;
Screen('TextStyle', window, 1);
Screen('TextSize', window, p.fontSize);
bbox = Screen('TextBounds', window, 'Loading Stimuli ... ');
newRect = CenterRectOnPoint(bbox, xCenter, yCenter);

Screen('DrawText', window, 'Loading Stimuli ... ', newRect(1), newRect(2), p.white);
Screen('Flip', window);

% Compute and store the center of the screen
p.xCenter = (p.sRect(3)-p.sRect(1))/2;
p.yCenter = (p.sRect(4)-p.sRect(2))/2;

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
    
    sinusoid = cos(2*pi*p.SF/w.ppd*(Xc.*cos(p.ori*(pi/180))+Yc.*sin(p.ori*(pi/180)))-p.phaseSteps(s));
      
    tmp = nan(size(Xc,1),size(Xc,2),2);
    tmp(:,:,1) = sinusoid.*p.maxContrastRGB+p.meanLum(1);
    tmp(:,:,2) = probe_aperture*255;
    probeText{s} = Screen('MakeTexture',window,tmp);
       
    tmp = nan(size(Xc,1),size(Xc,2),2);
    tmp(:,:,1) = sinusoid.*p.maxContrastRGB+p.meanLum(1);
    tmp(:,:,2) = upper_aperture*255;
    upperText{s} = Screen('MakeTexture',window,tmp);
    
    tmp = nan(size(Xc,1),size(Xc,2),2);
    tmp(:,:,1) = sinusoid.*p.maxContrastRGB+p.meanLum(1);
    tmp(:,:,2) = lower_aperture*255;
    lowerText{s} = Screen('MakeTexture',window,tmp);
    
    
end

toc

%% Determine stimulus sequence

stim.frameTimes = round(0:p.frameDur:p.runDur-p.frameDur,2)'; 
stim.nFrames = length(stim.frameTimes);

%% fixation targets

% determine when fixation targets occur
[stim.targTimes stim.nTargs] = createFixTargetSequence(p.runDur,p.frameDur,p.minTargDelay,p.maxTargDelay);

% create fixation target sequence
stim.targOnset = zeros(size(stim.frameTimes));
for t = 1:stim.nTargs
    [m,idx] = min(abs(stim.frameTimes-stim.targTimes(t))); % index closest frame
    stim.targOnset(idx) = 1; 
end
 
% determine when targets appear (3 frames = 250 ms)
stim.targ = zeros(size(stim.frameTimes));
for f = 1:stim.nFrames
     if f > 2 
         if stim.targOnset(f) == 1 || stim.targOnset(f-1) == 1 ||  stim.targOnset(f-2) == 1 % hard-coded
             stim.targ(f) = 1;
         end
     end 
end

% preallocate vectors to save response information
stim.resp = zeros(1,stim.nTargs);
stim.acc = zeros(1,stim.nTargs);
stim.rt = nan(1,stim.nTargs); 

%% stimulus sequence

% index when each stimulus is presented
stim.framesPerBlock = p.blockDur/p.frameDur;
stim.stimOn = []; 
for b = 1:p.blocksPerStim
    tmp = [zeros(1,stim.framesPerBlock),ones(1,stim.framesPerBlock)];
    stim.stimOn = [stim.stimOn tmp]; 
end
stim.stimOn = [stim.stimOn, zeros(1,stim.framesPerBlock)]; 

% mark the frames when each stimulus period starts
stim.stimStart = zeros(size(stim.stimOn)); 
stim.stimStart(1) = stim.stimOn(1);
for f = 2:stim.nFrames
   if stim.stimOn(f) == 1 & stim.stimOn(f-1) == 0
       stim.stimStart(f) = stim.stimOn(f);       
   end
end

%% Eye link setup
if p.eyeTrack == 1
    
    if p.runType == 0 % demo run
        edfFileFinal = [num2str(p.subID),'_FBA_CRF_WEDGES_Localizer_Demo',num2str(p.runNum),'.edf'];
    end
    if p.runType == 1 % regular run
        edfFileFinal = [num2str(p.subID),'_FBA_CRF_WEDGES_Localizer_Run',num2str(p.runNum),'.edf'];
    end
    
    edf_filename = 'tmp.edf';
    
    [el] = eyeTrackingOn_500Hz(window, edf_filename, p.bgRect, p.eyeRect, w.ppd, p.bgCol, p.txtCol);
    Screen('Flip', window);
end

%% wait for trigger
Screen('FillRect',window,p.bgCol,p.bgRect);
Screen('FillOval',window,p.white,p.fixRect); 
DrawFormattedText(window, '~', w.ppd, w.ppd, p.white, p.bgCol);
Screen('Flip', window);

% Create a queue which records keypresses from button box only
triggerKey = KbName('5%');
keylist = zeros(1,256); % keys for KbQueueCreate
keylist(p.keyPressNumbers) = 1; % only monitor p.keyPressNumbers

fprintf('Waiting for trigger... \n')

% Wait for scanner trigger
if p.scan == 1
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
presTime.frameTimes = presTime.runStart:p.frameDur:presTime.runEnd-p.frameDur;

% variables for recording psychtoolbox timing
time.vblstamp = nan(p.nFrames,1);
time.onset = nan(p.nFrames,1);
time.flipstamp = nan(p.nFrames,1);
time.missed = nan(p.nFrames,1);
time.runEnd = nan;

% initialize variables
targCnt = 1;
responseDeadline = nan;

%% begin task

for f = 1:p.nFrames
        
    %% determine phase of gratings
    upperPhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
    lowerPhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
    probePhase=round((rand(1,1)*(p.numPhaseSteps-1))+1);
        
    %% background
    Screen('FillRect',window,p.bgCol,p.bgRect);
    
    %% gratings
    if stim.stimOn(f) == 1      
        Screen('DrawTexture',window, upperText{upperPhase},[],p.attRect,0);
        Screen('DrawTexture',window, lowerText{lowerPhase},[],p.attRect,0);
        Screen('DrawTexture',window, probeText{probePhase},[],p.probeRect,0);
    end

     %% aperture and fixation
     if stim.targ(f) == 0
         Screen('FillOval',window,p.white,p.fixRect);
     else
         Screen('FillOval',window,p.fixTargCol,p.fixRect);
     end
       
    %% flip display
    
    Screen('DrawingFinished',window);
    
    % flip display
    [time.vblstamp(f), time.onset(f), time.flipstamp(f), time.missed(f)] = Screen('Flip',window,presTime.frameTimes(f));
    
    % send event marker
    if p.eyeTrack
        if stim.stimStart(f) == 0
            Eyelink('message','StimOff'); 
        end
        if stim.stimStart(f) == 1
           Eyelink('message','StimOn');
        end
    end
    
     if stim.targOnset(f)
         fprintf('Target %.0f, ',targCnt);
         PsychHID('KbQueueFlush',deviceNumber); % flush any response made before stimulus started
         targetStartTime = GetSecs;
         responseDeadline = targetStartTime + p.responseWindow;
     end
        
    %% monitor for response

    while 1
        
        % break loop after frame duration
        if GetSecs >= presTime.frameTimes(f)+p.frameDur
            break
        end
        
        % if responseDeadline exceeded....
        if GetSecs > responseDeadline
            fprintf('Acc = %.0f\n',stim.acc(targCnt));
            responseDeadline = nan; % set responseDeadline to a nan (so targCnt won't advnace until next trial
            targCnt = targCnt + 1;  % advance target counter
        end
        
        % check for response if haven't hit the response deadline
        if GetSecs <= responseDeadline
            
            [pressed, keyIndex] = PsychHID('KbQueueCheck', deviceNumber);
            
            if pressed
                
                rt = GetSecs - targetStartTime; % calculate rt
                
                [m, keyCode] = max(keyIndex);
                                
                PsychHID('KbQueueFlush',deviceNumber);
                
                % save response time
                stim.rt(targCnt) = rt;
                
                % save response
                stim.resp(targCnt) = 1;
                
                % mark as response made
                stim.acc(targCnt) = 1;
                
            end
            
        end       
        
    end
    
    
end

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
    DrawFormattedText(window, Text, 'center', yCenter-w.ppd*2, p.white);
    for n = 1:numel(targLocations)
        textRect = CenterRectOnPoint(Screen('TextBounds', window, targText{n}), xCenter-targLocations(n), yCenter);
        DrawFormattedText(window, targText{n}, textRect(1), yCenter-w.ppd*1, p.white);
        Screen('FillOval', window, p.white, CenterRectOnPoint([0 0 round(w.ppd) round(w.ppd)], xCenter+targLocations(n), yCenter+w.ppd/2));
        Screen('FillOval', window, p.bgCol, CenterRectOnPoint([0 0 round(w.ppd)-4 round(w.ppd)-4], xCenter+targLocations(n), yCenter+w.ppd/2));
        numRect = CenterRectOnPoint(Screen('TextBounds', window, num2str(n)), xCenter+targLocations(n), yCenter+w.ppd/2);
         DrawFormattedText(window, num2str(n), numRect(1), numRect(2), p.white);
    end
    Screen('FillOval', window, [255 0 0], CenterRectOnPoint([0 0 round(w.ppd) round(w.ppd)], xCenter+targLocations(startloc), yCenter+w.ppd/2));
    Screen('FillOval', window, p.bgCol, CenterRectOnPoint([0 0 round(w.ppd)-4 round(w.ppd)-4], xCenter+targLocations(startloc), yCenter+w.ppd/2));
    numRect = CenterRectOnPoint(Screen('TextBounds', window, num2str(startloc)), xCenter+targLocations(startloc), yCenter+w.ppd/2);
    DrawFormattedText(window, num2str(startloc), numRect(1), numRect(2), p.white);
    Screen('Flip', window);
end

wake.WakefullnessRating = startloc;
wake.WakefullnessLabels = fliplr(targText);

%% Clean up

% save the data
save(fName,'p','w','stim','presTime','time','wake');

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
% if p.scan == 1
%     Datapixx('RegWrRd');
%     Datapixx('close');
% end

PsychHID('KbQueueStop', deviceNumber);
ShowCursor;
Priority(0);

% print performance

accuracy = 100*length(stim.acc(stim.acc == 1))/stim.nTargs; 
hits = length(stim.acc(stim.acc == 1)); 
fprintf('Percent Correct = %.1f, Hits = %.0f out of %.0f \n',accuracy,hits,stim.nTargs)