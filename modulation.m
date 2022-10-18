function [x1, x2] = modulation(data, index, numSymbols)

            i = sqrt(-1);
            
            bits = char(string(data(index) + string(data(index + 1))));
            n = bin2dec(bits);
            x1 = exp(i*pi*((2*n+1)/numSymbols));
            
            bits = char(string(data(index + 2) + string(data(index + 3))));
            n = bin2dec(bits);
            x2 = exp(i*pi*((2*n+1)/numSymbols));
           
end

