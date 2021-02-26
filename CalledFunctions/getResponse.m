function resp = getResponse(p,w,cueCol,cueValidity)

respKey = KbName('SPACE'); % response starts as 32

Screen('FillRect',w,p.bgCol,p.bgRect);
% Screen('DrawTexture',w,postCue,[],p.cueRect);
Screen('FillOval',w,cueCol,p.fixRect);
text = 'Target?';
tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-p.questionDist];
Screen('DrawText', w, text, tCenter1(1), tCenter1(2), p.txtCol);
if cueValidity == 0
    text = 'Invalid cue';
    tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter+p.questionDist/2];
    Screen('DrawText', w, text, tCenter1(1), tCenter1(2), 0);
end
Screen('DrawingFinished',w);
Screen('Flip',w);

while KbCheck; end % Wait until all keys are released before checking for input

while 1
    [ keyIsDown, seconds, keyCode ] = KbCheck; % Check the state of the keyboard
    if keyIsDown
        kp = find(keyCode); % get the key code
        if ismember(kp, p.acceptKeys) % only proceed if acceptable key pressed. Otherwise, begin loop again.
            
            while KbCheck; end % wait for keys to be released until moving on
            
            % if more than one key has been recorded randomly pick one.
            if length(kp) > 1
                kp = randsample(kp,1);
            end
            
            if respKey ~= KbName('SPACE') && kp == KbName('SPACE') % if a non-spacebar response is recorded and spacebar has been pressed, exit function
                break
            else
                respKey = kp; % else update the response.
            end
            
            % Show the response (only if it's not spacebar)
            if respKey ~= KbName('SPACE')
                
                if respKey == p.keys.yes
                    respText = 'Yes';
                    resp = 1;
                end
                if respKey == p.keys.no
                    respText = 'No';
                    resp = 0;
                end
                Screen('FillRect',w,p.bgCol,p.bgRect);
                Screen('FillOval',w,cueCol,p.fixRect);
                text = 'Target?';
                tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-p.questionDist];
                Screen('DrawText', w, text, tCenter1(1), tCenter1(2), p.txtCol);
                tCenter2 = [p.xCenter-RectWidth(Screen('TextBounds', w, respText))/2 p.yCenter-p.answerDist];
                Screen('DrawText', w, respText, tCenter2(1), tCenter2(2), p.txtCol);
                if cueValidity == 0
                    text = 'Invalid cue';
                    tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter+p.questionDist/2];
                    Screen('DrawText', w, text, tCenter1(1), tCenter1(2), 0);
                end
                Screen('DrawingFinished',w);
                Screen('Flip',w);
            end
            
        else
            continue
        end
    end
end


