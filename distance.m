clear;clc; close all; clear functions;
addpath(genpath(pwd));

distName = "hamming";

%delete(gcp('nocreate'));
%parpool("Processes", 8);

outputFolder = "Data/"+distName+"/distChr/";
createFolder(outputFolder);


for chr=1:1:19
    % Computing gene distance per chromosome chr
    %distGene(chr, distName);

    % Chromosome distance folder
    outputMatrixFolder = "Data/"+distName+"/chr"+chr+"/distances";
    
    % Computing chromosome distance
    chromosomeDistance = meanMatrix(outputMatrixFolder);
    
    % Saving chromosome distance
    writematrix(chromosomeDistance, outputFolder+"chr"+chr+".txt");

    disp("Chromosome completed");

end

globalDistance = meanMatrix(outputFolder);
writematrix(globalDistance, "Data/"+distName+"/globalDist.txt");



