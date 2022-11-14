function [x1, x2] = modulation(data, index, numSymbols)

            i = sqrt(-1);
            symbolSize = log2(numSymbols);
            
            bitsVec = data(index:index + symbolSize - 1);
            bits = "";
            for j = 1:symbolSize
                bits = bits + string(bitsVec(j));
            end
            bits = char(bits);
            n = bin2dec(bits);
            x1 = exp(i*pi*((2*n+1)/numSymbols));
            
            bitsVec = data(index + symbolSize:index + symbolSize* 2 - 1);
            bits = "";
            for j = 1:symbolSize
                bits = bits + string(bitsVec(j));
            end
            bits = char(bits);
            n = bin2dec(bits);
            x2 = exp(i*pi*((2*n+1)/numSymbols));
           
end

