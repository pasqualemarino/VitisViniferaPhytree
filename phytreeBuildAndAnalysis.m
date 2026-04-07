clear;clc; close all; clear functions;

addpath(genpath(pwd));

%% Tree build
function TreeBuildNJ()
    names = string(readcell("Data/names.txt", 'NumHeaderLines', 0,  'TextType', 'string', "Delimiter",""));
    distGlobal = readmatrix("Data/hamming/DistanzaTraIndividui.txt");
    GlobalTree = seqneighjoin(distGlobal, 'equivar', names);
    phytreewrite("Data/hamming/treeNJ/global.tree", GlobalTree);
    for i=1:19
        matrixChr = readmatrix("Data/hamming/distChr/chr"+i+".txt");
        treeChr = seqneighjoin(matrixChr, 'equivar', names);
        phytreewrite("Data/hamming/treeNJ/Chr"+i+".tree", treeChr);
    end
end

function TreeBuildME()
    names = string(readcell("Data/names.txt", 'NumHeaderLines', 0,  'TextType', 'string', "Delimiter",""));
    distGlobal = readmatrix("Data/hamming/DistanzaTraIndividui.txt");
    GlobalTree = METree(distGlobal);
    phytreewrite("Data/hamming/treeME/global.tree", GlobalTree);
    for i=1:19
        matrixChr = readmatrix("Data/hamming/distChr/chr"+i+".txt");
        treeChr = METree(matrixChr);
        phytreewrite("Data/hamming/treeME/Chr"+i+".tree", treeChr);
    end
end

function Tree = METree(distMatrix)
    matrix = squareform(distMatrix, "tomatrix");
    writematrix(matrix, "R/matrixME.txt");
    system('/usr/local/bin/R CMD BATCH R/ME.R R/outME.txt');
    Tree = phytreeread("R/TreeME.tree");
end

%% Topology NJ: Robinson-Fould distance
function distRF = robinsonFouldDistanceMatrix()
    distRF = zeros(20, 20);
    GlobalTree = phytreeread("Data/hamming/treeNJ/global.tree");
    for i = 1:19
        Tree1 = phytreeread("Data/hamming/treeNJ/Chr"+i+".tree");
        distRF(20, i) = robinsonFoulds(Tree1, GlobalTree);
        distRF(i, 20) = distRF(20, i);
        for j = i+1:19
            Tree2 = phytreeread("Data/hamming/treeNJ/Chr"+j+".tree");

            distRF(j, i) = robinsonFoulds(Tree1, Tree2);
            distRF(i, j) = distRF(j, i);
        end
    end
end


% Computing RB distance of a pair of tree
function n=robinsonFoulds(tree1, tree2)
    phytreewrite("R/tree1.tree", tree1);
    phytreewrite("R/tree2.tree", tree2);
    system('/usr/local/bin/R CMD BATCH R/rf.R R/outRF.txt');
    n = readmatrix("R/dist.txt");
end

%% Patristic distance NJ
function patristicDist = patristicDistanceMatrix()
    patristicDist = zeros(20, 20);
    GlobalTree = phytreeread("Data/hamming/treeNJ/global.tree");
    distGTree = pdist(GlobalTree, "squareform", true);
    nameGlobal = get(GlobalTree, "LeafNames");
    for i=1:19
        treeChr1 = phytreeread("Data/hamming/treeNJ/Chr"+i+".tree");
        distChr1 = pdist(treeChr1, "squareform", true);
        nameChr1 = get(treeChr1, "LeafNames");
        [~, orderCorrect] = ismember(nameChr1, nameGlobal);
        distGTree = distGTree(orderCorrect, orderCorrect);
        patristicDist(20, i) = norm(distChr1-distGTree, "fro");
        patristicDist(i, 20) = patristicDist(20, i);
        for j=i+1:19
            treeChr2 = phytreeread("Data/hamming/treeNJ/Chr"+j+".tree");
            distChr2 = pdist(treeChr2, "squareform", true);
            nameChr2 = get(treeChr2, "LeafNames");
            [~, orderCorrect] = ismember(nameChr1, nameChr2);
            distChr2 = distChr2(orderCorrect, orderCorrect);
            patristicDist(j, i) = norm(distChr1 - distChr2, "fro");
            patristicDist(i, j) = patristicDist(j, i);
        end
    end
end

%% Topology NJ vs ME: Robinson-Fould distance
function dist = distRF_NJvsME()
    dist = zeros(20, 1);
    for i=1:19
        TreeNJ = phytreeread("Data/hamming/treeNJ/Chr"+i+".tree");
        TreeME = phytreeread("Data/hamming/treeME/Chr"+i+".tree");
        dist(i) = robinsonFoulds(TreeNJ, TreeME);
    end
    TreeNJ = phytreeread("Data/hamming/treeNJ/global.tree");
    TreeME = phytreeread("Data/hamming/treeME/global.tree");
    dist(20) = robinsonFoulds(TreeNJ, TreeME);
end

%% Patristic distance NJ vs ME
function dist = patristicDist_NJvsME() %Sistemare l'ordine della seconda matrice
    dist = zeros(20, 1);
    for i=1:19
        TreeNJ = phytreeread("Data/hamming/treeNJ/Chr"+i+".tree");
        TreeME = phytreeread("Data/hamming/treeME/Chr"+i+".tree");
        distNJ = pdist(TreeNJ, "squareform", true);
        distME = pdist(TreeME, "squareform", true);
        NJNames = get(TreeNJ, "LeafNames");
        MENames = get(TreeME, "LeafNames");
        [~, orderCorrect] = ismember(NJNames, MENames);
        distME = distME(orderCorrect, orderCorrect);
        diff = distNJ-distME;
        dist(i) = norm(diff, "fro");
    end
    TreeNJ = phytreeread("Data/hamming/treeNJ/global.tree");
    TreeME = phytreeread("Data/hamming/treeME/global.tree");
    distNJ = pdist(TreeNJ, "squareform", true);
    distME = pdist(TreeME, "squareform", true);
    NJNames = get(TreeNJ, "LeafNames");
    MENames = get(TreeME, "LeafNames");
    [~, orderCorrect] = ismember(NJNames, MENames);
    distME = distME(orderCorrect, orderCorrect);
    diff = distNJ-distME;
    dist(20) = norm(diff, "fro");
end
