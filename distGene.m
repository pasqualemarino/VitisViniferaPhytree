function distGene(chr, distanceName)

    outputFolder = "Data/"+distanceName+"/chr"+chr+"/";
    matrixOutputFolder = "Data/"+distanceName+"/chr"+chr+"/matrici_distanze/";
    
    % Folder checking
    createFolder("Data/"+distanceName+"/chr"+chr);
    createFolder(matrixOutputFolder);
    
    % All possible gene in chromosome chr
    GeneFolder = "../Dati/all_allele_sequences/chr"+chr+"/";
    genes = dir(fullfile(GeneFolder, '*.txt'));
    numberGenes = length(genes);
    
    for i=1:numberGenes
    
        % Retrieving gene's name
        [~, geneName, ~] = fileparts(genes(i).name);
        %disp("Analizzo gene "+nomeGene + " -  chr "+chr+" : gene "+i+" su " + numeroGeni);
    
        % haplotype folder
        hapFolder = "../Dati/all_gemes_nofilter-2/ch"+chr+"/"+geneName+"_serie_allelica.txt";
    
        % Retrieving alleles
        Alleles = readmatrix(GeneFolder+genes(i).name,  'OutputType', 'char');
        Alleles = Alleles(2:end, :);
        
        % Check on alleles
        if ~any(ismissing(Alleles), "all")
            % Computing alleles distance
            distAll = distAlleles(Alleles, distanceName);
            
            % Retrieving haplotypes
            haplotype = readcell(hapFolder, 'NumHeaderLines', 0,  'TextType', 'string');
            haplotype = string(haplotype);
        
            % Removing header and indefity 81
            haplotype = haplotype(2:end, 1:3);
            haplotype( 81, :) = [];
            
            % Computing gene distance
            obj = distSingleGene();
            obj.distAllelesMatrix = distAll;
            
            try
                DistGene = obj.distGeneMatrix(haplotype);
                distVector = squareform(DistGene, "tovector");
                writematrix(distVector, matrixOutputFolder+geneName+".txt");
            catch exception
                disp("Error in computing gene distance "+geneName+" ("+i+")");
            end
        end
    end
end