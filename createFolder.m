function createFolder(folder)
    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end