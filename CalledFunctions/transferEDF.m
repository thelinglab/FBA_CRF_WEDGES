function transferEDF(edf_filename,edf_filename_final)
% this function copies and renames file from eye tracking host to stim pres
% machine.

try
    fprintf('Receiving data file ''%s''\n', edf_filename);
    status=Eyelink('ReceiveFile');
    if status > 0
        fprintf('ReceiveFile status %d\n', status);
    end
    if 2==exist(edf_filename, 'file')      
        
        inputFullFileName = fullfile(pwd, edf_filename);
        outputFullFileName = fullfile(pwd, edf_filename_final);
        copyfile(inputFullFileName, outputFullFileName);
        
    end
catch
    fprintf('Problem receiving data file ''%s''\n', edf_filename);
end