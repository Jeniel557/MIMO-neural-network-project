close all;
clear all;
clc;

i = sqrt(-1);
SNR = 6;           % Signal-to-noise ratio (dB)
rng default

hAWGN = comm.AWGNChannel(...
    'NoiseMethod','Signal to noise ratio (SNR)',...
    'SNR',SNR,...
    'SignalPower',1);

%% Transmission

numBits = 4000;  
data_train = randi([0 1],numBits,1);
data_test = randi([0 1],numBits,1);
numSymbols = 4;
symbolSize = log2(numSymbols);



for phase = 1:2
%     phase
    if phase == 1 % We are generating training data
        data = data_train;
    else
        data = data_test;
    end
      
    modData = [];
    recivedSymbols = [];
    recived = [];
    prediction = [];
        
    for j = 1:2*symbolSize:length(data)
    
        % Modulate the data into symbols.
            bits = char(string(data(j) + string(data(j + 1))));
            n = bin2dec(bits);
            x1 = exp(i*pi*((2*n+1)/numSymbols));
            
            bits = char(string(data(j + 2) + string(data(j + 3))));
            n = bin2dec(bits);
            x2 = exp(i*pi*((2*n+1)/numSymbols));
            
            modData = [modData; [x1; x2]];
            modData_real = [real(modData); imag(modData)];

        % Encode the symbols with the Alamoteee
            t1 = [x1; x2];
            t2 = [conj(-x2); conj(x1)];
            
        % Transmission of the symbols

            % channel matrix 
            scale = 0.75; % set the scale to use as well as the size of the H matrix.
            h_rows = 2; 
            h_columns = 2;
            H = KnownChannelMatrixCreation(scale,h_rows,h_columns) ;
            
            R = [H*t1; conj(H*t2)];
          
            % noise
            R = step(hAWGN,R);
            
       % Channel estimation
        
        H_est = channelEstimation(R,t1,t2);

       % Decoder     

            H_bar = [conj(H_est(1,2)), -conj(H_est(1,1)); conj(H_est(2,2)), -conj(H_est(2,1))];
            decoding_H = [H_est; H_bar];  
            x_decoded = (inv(decoding_H' * decoding_H)* decoding_H')*R;
            x_decoded(2) = conj(x_decoded(2));
             
            recivedSymbols = [recivedSymbols; x_decoded];
            x_estimate_real = [real(recivedSymbols); imag(recivedSymbols)];
            
            if phase == 2
                
                x_decoded(2) = conj(x_decoded(2));
                x_estimate_real = [real(x_decoded); imag(x_decoded)];
                estimate = predict(net,x_estimate_real);
                
                estimate_symbols = estimate(1) + estimate(3)*i;
                prediction = [prediction; estimate_symbols];
                
                estimate_symbols = estimate(2) + estimate(4)*i;
                prediction = [prediction; estimate_symbols];
            end
    end
    
    if phase == 1
        net = fitcnet(x_estimate_real,modData_real,"LayerSizes",[10 10 10])
        figure;
        plot(net.TrainingHistory.TrainingLoss)
    end  
end

prediction - modData;


%% convert back to bits

        % figure out which symbol maps to which bits.
       
       mapping =  mapSymbolsToBits(modData, data);



        % convert the symbols back into bits.
        output_str = "";
        for i = 1:size(prediction,1)
            output_str = output_str + mapping(find(mapping == string(prediction(i))) + 1);
        end
        output_str = char(output_str);


        % convert from a string back into a vector.
        outputVector = [];
        for i = 1:ceil(size(data,1))
            outputVector = [outputVector; str2num(output_str(i))];
        end

        numberOfBitErrors = sum(abs(outputVector -  data))

        BER = numberOfBitErrors/size(data,1)
        
        
        
        
      %%  Saving data to a file
%      
% Create two columns of data
A  = rand(10,1);
B = rand(10,1);
% Create a table with the data and variable names
T = table(A, B, 'VariableNames', { 'A', 'B'} )
% Write data to text file
writetable(T, 'MyFile.txt')
        
        
%% Graphs

figure;
plot(modData,'xr','MarkerSize',15)
hold on
plot(recivedSymbols, 'ob')
xlim([-2 2])
ylim([-2 2])
title("decoded vs transmitted symbols")

figure;
plot(modData,'xr','MarkerSize',15)
hold on
plot(prediction, 'ob')
xlim([-2 2])
ylim([-2 2])
title("decoded vs transmitted symbols")