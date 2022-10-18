function [outputVector] = convertToBits(modData, data, prediction )
%MAPSYMBOLSTOBITS Summary of this function goes here
%   Detailed explanation goes here

        mapping = [];
        possibleSymbols = unique(modData , 'first');
        
%         mapSymbolsToBits

        for i = 1:size(possibleSymbols,1)
            index = 1;
            while index < size(modData,1)

                if(modData(index) == possibleSymbols(i))
                    start = (index -1)*2 + 1;
                    bits = string(data(start)) + string(data(start + 1));
                    mapping = [mapping,modData(index)];
                    mapping = [mapping,bits];
                    index = size(modData,1);
                end
                index = index + 1;
            end
        end
        
        % convert the symbols back into bits.
        output_str = "";
        for i = 1:size(prediction,1)
            output_str = output_str + mapping(find(mapping == string(prediction(i))) + 1);
        end
        output_str = char(output_str);
        outputVector = [];
        for i = 1:ceil(size(data,1))
            outputVector = [outputVector; str2num(output_str(i))];
        end


end

