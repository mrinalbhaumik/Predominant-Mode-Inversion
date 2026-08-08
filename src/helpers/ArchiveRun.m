function archivePath = ArchiveRun(mainScriptName, targetFolder)
% Archives every .m file a script depends on into a timestamped folder, creating a
% reproducible snapshot of the code used for a run.
% Author      : Mrinal Bhaumik
% Affiliation : Utah State University, 2026

    % Create timestamped subfolder
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    archivePath = fullfile(targetFolder, ['archive_', timestamp]);
    mkdir(archivePath);

    % Get list of all dependent .m files
    usedFiles = matlab.codetools.requiredFilesAndProducts(mainScriptName);

    % Copy them to archive folder
    for i = 1:length(usedFiles)
        srcf = usedFiles{i};
        relPath = strrep(srcf, pwd, '');
        dest = fullfile(archivePath, relPath);
        destFolder = fileparts(dest);
        if ~exist(destFolder, 'dir')
            mkdir(destFolder);
        end
        copyfile(srcf, dest);
    end

    fprintf('Script and dependencies archived to: %s\n', archivePath);
end
