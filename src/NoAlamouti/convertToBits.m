function [output] = convertToBits(prediction,numSymbols)
%MAPSYMBOLSTOBITS Summary of this function goes here
%   Detailed explanation goes here

        i = sqrt(-1);
        symbolSize = log2(numSymbols);
        % convert the symbols back into bits.
           
           outputVector = [];
           for j = 1:length(prediction)
               
                n = (numSymbols*i*log(prediction(j)))/(-2*pi) - 0.5;
                if n < 0
                    n = n + numSymbols;
                end
                n = round(n,1);
                
                if n == numSymbols
                    n = 0;
                end
                outputVector = [outputVector;real(n)];
           end

           output_str = "";
           for j = 1:length(outputVector)

               temp = string(dec2bin(outputVector(j)));
               while length(char(temp)) < symbolSize
                   temp = "0" + temp;
               end
               output_str = output_str  + temp;
           end
       
        output_str = char(output_str);
        output = [];
        for j = 1:length(output_str)
            output =[output;str2num(output_str(j))];
        end
end