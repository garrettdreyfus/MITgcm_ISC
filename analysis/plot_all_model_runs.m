rootdir = '/home/garrett/Projects/MITgcm_ISC/experiments/'
foldernames =  GetSubDirsFirstLevelOnly(rootdir)
for i = 1:numel(foldernames)
    experimentnames =  GetSubDirsFirstLevelOnly([rootdir foldernames{i}]);
    foldernames{i}
    for j = 1:numel(experimentnames)
        try
            plot_model([rootdir foldernames{i}], experimentnames{j},['geometries/' foldernames{i} '|' experimentnames{j}],false);
        catch
            [rootdir foldernames{i}]
        end
    end
end


function [subDirsNames] = GetSubDirsFirstLevelOnly(parentDir)
    % Get a list of all files and folders in this folder.
    files = dir(parentDir);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subDirs = files(dirFlags);
    subDirsNames = cell(1, numel(subDirs) - 2);
    for i=3:numel(subDirs)
        subDirsNames{i-2} = subDirs(i).name;
    end
end
