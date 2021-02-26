function wait_for_spacebar(p,w)
% wait for spacebar press

while 1
    [keyIsDown,secs,keyCode]=KbCheck;
    if keyIsDown
        kp = find(keyCode);
        if kp == p.keys.space  % if spacebar
            break;
        end
        if kp == p.keys.esc % if escape
            LoadIdentityClut(w); % restore identity clut
            Screen('CloseAll'); % close screen
            return
        end
    end
end
WaitSecs(.1);