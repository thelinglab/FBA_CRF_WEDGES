%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      EYE LINK SETUP                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calibration_area_proportion - if I use this, don't have to specify exact
% location of calibration/validation points, instead it scales the screen.

% download edfmex from eyelink forum to import eyelink output directly into
% matlab

function [el] = eyeTrackingOn(w, edf_filename, rect, eyeRect, pixelPerDeg, bgCol, txtCol)

% Provide Eyelink with details about the graphics environment
% and perform some initializations. The information is returned
% in a structure that also contains useful defaults
% and control codes (e.g. tracker state bit and Eyelink key values).
el=EyelinkInitDefaults(w);%%returns values
% setup the proper calibration foreground and background colors
el.backgroundcolour = bgCol;
el.msgfontcolour = txtCol;
el.foregroundcolour = BlackIndex(el.window);
el.calibrationtargetcolour= WhiteIndex(el.window);
%initialize buffers to background color
Screen('FillRect',w, bgCol);
Screen('Flip', w,0,1);

Screen('FillRect',w, bgCol);
Screen('Flip', w,0,1);

EyelinkUpdateDefaults(el);

% STEP 4
% Initialization of the connection with the Eyelink Gazetracker.
% exit program if this fails.
if ~EyelinkInit(0)
    fprintf('Eyelink Init aborted.\n')
    % Shutdown Eyelink:
    Eyelink('Shutdown');
    return;
end

connected=Eyelink('IsConnected');%Just to verify connection

[v, vs]=Eyelink('GetTrackerVersion');

%edf_filename =['Pupil_', [num2str(SUBJECT) '_' num2str(runN) '_'num2str(Session)], '.edf'];
%edf_filename =['P', [num2str(SUBJECT) num2str(runN) num2str(Session)], '.edf'];
% edf_filename =[num2str(observer) '.edf'];
% [status] = Eyelink('Openfile', edf_filename);
fprintf('Running experiment on a ''%s'' tracker.\n', vs );

% open file to record data to
tempeye = Eyelink('Openfile', edf_filename);%%doesn't return
%EyelinkUpdateDefaults(el);
if tempeye~=0
    fprintf('Cannot create EDF file ''%s'' ', edf_filename);
    Eyelink( 'Shutdown');
    return;
end

% CHANGE HOST PC PARAMETERS HERE
% SET UP TRACKER CONFIGURATION
% Setting the proper recording resolution, proper calibration type,
% as well as the data file content;
Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, (rect(3)-1), (rect(4)-1)); % scr_r(3) = swidth; scr_r(4) = sheight
Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, (rect(3)-1), (rect(4)-1));

% set calibration type.
Eyelink('command', 'calibration_type = HV5');

% set pupil Tracking model in camera setup screen no = centroid. yes = ellipse
Eyelink('command', 'use_ellipse_fitter = yes');

% Make sure we're recording the gaze
Eyelink('command', 'recording_parse_type = GAZE');
Eyelink('command', 'sample_rate = 500');

%%%%%%%%%%%%%%%%%%%custom calibration

% you must send this command with value NO for custom calibration
% you must also reset it to YES for subsequent experiments
Eyelink('command', 'generate_default_targets = NO');

% STEP 5.1 modify calibration and validation target locations
SCREEN_NO = max(Screen('Screens'));
[swidth, sheight]=Screen('WindowSize', SCREEN_NO);
height=sheight;
width=swidth;
caloffset=round(7*pixelPerDeg);%%%%% 5 point

Eyelink('command','calibration_samples = 6');
Eyelink('command','calibration_sequence = 0,1,2,3,4,5');

Eyelink('command','calibration_targets = %d,%d %d,%d %d,%d %d,%d %d,%d',...
round(width/2),round(height/2),  round(width/2),round(height/2-caloffset),  round(width/2),round(height/2 + caloffset),  round(width/2 -caloffset),round(height/2),  round(width/2 +caloffset),round(height/2) );

Eyelink('command','validation_samples = 6');
Eyelink('command','validation_sequence = 0,1,2,3,4,5');
Eyelink('command','validation_targets = %d,%d %d,%d %d,%d %d,%d %d,%d',...
round(width/2),round(height/2),  round(width/2),round(height/2-caloffset),  round(width/2),round(height/2 + caloffset),  round(width/2 -caloffset),round(height/2),  round(width/2 +caloffset),round(height/2));

%%% 9 point calibration
% caloffset=round(7.5*pixelPerDeg);
% Eyelink('command','calibration_samples = 10');
% Eyelink('command','calibration_sequence = 0,1,2,3,4,5,6,7,8,9');
% Eyelink('command','calibration_targets = %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d',...
% round(swidth/2),round(sheight/2),  round(swidth/2),round(sheight/2)-caloffset,  round(swidth/2),round(sheight/2) + caloffset,  round(swidth/2) -caloffset,round(sheight/2),  round(swidth/2) +caloffset,round(sheight/2),...
% round(swidth/2)-caloffset, round(sheight/2)- caloffset, round(swidth/2)-caloffset, round(sheight/2)+ caloffset, round(swidth/2)+caloffset, round(sheight/2)- caloffset, round(swidth/2)+caloffset, round(sheight/2)+ caloffset);
% Eyelink('command','validation_samples = 9');
% Eyelink('command','validation_sequence = 0,1,2,3,4,5,6,7,8,9');
% Eyelink('command','validation_targets = %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d %d,%d',...
% round(swidth/2),round(sheight/2),  round(swidth/2),round(sheight/2)-caloffset,  round(swidth/2),round(sheight/2) + caloffset,  round(swidth/2) -caloffset,round(sheight/2),...
% round(swidth/2) +caloffset,round(sheight/2),...
% round(swidth/2)-caloffset, round(sheight/2)- caloffset, round(swidth/2)-caloffset, round(sheight/2)+ caloffset, round(swidth/2)+caloffset, round(sheight/2)- caloffset, round(swidth/2)+caloffset, round(sheight/2)+ caloffset);

Eyelink('command', 'select_parser_configuration = 0');

Eyelink('command', 'saccade_acceleration_threshold = 8000');

Eyelink('command', 'saccade_velocity_threshold = 30');
 
Eyelink('command', 'saccade_motion_threshold = 0.0');
 
Eyelink('command', 'saccade_pursuit_fixup = 60');
 
Eyelink('command', 'fixation_update_interval = 0');

% set EDF file contents
Eyelink('command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS');

% set link data (used for gaze cursor)
Eyelink('command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');

% make sure we're still connected.
if Eyelink('IsConnected')~=1
    Eyelink( 'Shutdown');
    return;
end;

% Hide the mouse cursor;
Screen('HideCursorHelper', w);

% Calibrate the eye tracker
EyelinkDoTrackerSetup(el);
eye_used = Eyelink('EyeAvailable');

% Must be offline to draw to EyeLink screen
Eyelink('Command', 'set_idle_mode');

% clear tracker display and draw box at fix point
Eyelink('Command', 'clear_screen 0')

Eyelink('command', 'draw_box %d %d %d %d 15',eyeRect(1),eyeRect(2),eyeRect(3),eyeRect(4));

Eyelink('Command', 'set_idle_mode');
WaitSecs(0.01);%.05