classdef distSingleGene
   properties
      distAllelesMatrix
   end
   methods
        
        function D = distGeneMatrix(obj, dataMatrix)
            numI = size(dataMatrix, 1);
        
            D = zeros(numI, numI);
            
            for i=1:1:numI
                I1_A = dataMatrix(i, 2);
                I1_B = dataMatrix(i, 3);
                for j=i+1:1:numI
                    I2_A = dataMatrix(j, 2);
                    I2_B = dataMatrix(j, 3);
                    D(j, i) = obj.distMultipleAlleles(dataMatrix(i, 2), dataMatrix(i, 3), dataMatrix(j, 2), dataMatrix(j, 3) );
                end

            
            end
        end

        function k = distMultipleAlleles(obj, I1_A, I1_B, I2_A, I2_B)
            if I1_A=="not_found" && I1_B~="not_found"
                I1_A = I1_B;
            end
            if I1_B=="not_found" && I1_A~="not_found"
                I1_B = I1_A;
            end
            if I2_A=="not_found" && I2_B~="not_found"
                I2_A = I2_B;
            end
            if I2_B=="not_found" && I2_A~="not_found"
                I2_B = I2_A;
            end
            
            % Retrieving all possibile alleles
            indexs = @(x) str2double( extractAfter(x, 4) );
            alleles1_A = arrayfun(indexs, split(I1_A, "|"));
            alleles1_B = arrayfun(indexs, split(I1_B, "|"));
            alleles2_A = arrayfun(indexs, split(I2_A, "|"));
            alleles2_B = arrayfun(indexs, split(I2_B, "|"));
           
            % Computing all possible genotype
            [V1, V2] = ndgrid(alleles1_A, alleles1_B);
            potentialI1 = [V1(:), V2(:)];
            potentialI1 = sort(potentialI1, 2);

            [V3, V4] = ndgrid(alleles2_A, alleles2_B);
            potentialI2 = [V3(:), V4(:)];
            potentialI2 = sort(potentialI2, 2);
            

            [genotype1, ~, i1] = unique(potentialI1, 'rows');
            weight1 = accumarray(i1, 1);

            [genotype2, ~, i2] = unique(potentialI2, "rows");
            weight2 = accumarray(i2, 1);
        
            % Computing distance
            sumA = 0;
            for i=1:1:size(genotype1, 1)
                for j=1:1:size(genotype2, 1)
                    dist = weight1(i) * weight2(j) * obj.distanceSingleAlleles(genotype1(i, 1), genotype1(i, 2), genotype2(j, 1), genotype2(j, 2));
                    sumA = sumA + dist;
                end
            end 
            numberSums = sum(weight1) * sum(weight2);
            k = sumA / numberSums; 
           
        end

        % Computing single allele distance
        function k = distanceSingleAlleles(obj, I1_A, I1_B, I2_A, I2_B)
            sum1 = [obj.distAllelesMatrix(I1_A, I2_A)] + [obj.distAllelesMatrix(I1_B, I2_B)];
            sum2 = [obj.distAllelesMatrix(I1_A, I2_B)] + [obj.distAllelesMatrix(I1_B, I2_A)];
            if sum1<sum2
                k = sum1;
            else
                k = sum2;
            end
        end
        
        function i= AlleleToIndex(obj, alleles)
            i = str2double( extractAfter(alleles, 4) );
        end
   end
end