function [sym] = argmax(x,numSymbols)

    i = sqrt(-1);
    
    min = 1000;
    index = 0;
    for j = 1:numSymbols
        sym = exp(i*pi*((2*j+1)/numSymbols));
        dist = sqrt((real(sym) - real(x))^2  + (imag(sym) - imag(x))^2);

        if dist < min
            min = dist;
            index = j;
        end
    end

    sym = exp(i*pi*((2*index+1)/numSymbols));
end

