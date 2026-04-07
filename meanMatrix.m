% Computing mean of matrix inside folder
function M = meanMatrix(folder)

    % Retrieving files
    files = dir(fullfile(folder, '*.txt'));
    fileCount = length(files);

    % Computing mean
    sum = 0;
    for i = 1:fileCount
        fileName = fullfile(folder, files(i).name);
        currentMatrix = readmatrix(fileName);
        sum = sum + currentMatrix;
        if any(isnan(currentMatrix), "all")
            disp(fileName);
        end
    end
    M = sum / fileCount;
end