close all;
clear all;
clc;

i = sqrt(-1);
SNR = 6;
rx = 2; 
tx = 2;


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
            [x1, x2] = modulation(data, j, numSymbols);
            
            modData = [modData; [x1; x2]];
            modData_real = [real(modData); imag(modData)];

        % Encode the symbols with the Alamoteee
            t1 = [x1; x2];
            t2 = [conj(-x2); conj(x1)];
            
        % Transmission of the symbols
        
            % channel matrix 
            scale = 0.75; % set the scale to use as well as the size of the H matrix.
            H = KnownChannelMatrixCreation(scale,rx,tx) ;
            
            % Recieved values.
            R = [H*t1; conj(H*t2)];
            R = step(hAWGN,R);
            
       % Channel estimation
            H_est = channelEstimation(R,t1,t2);

       % Decoder    
            x_decoded = decoder(H_est, R);
            
            recivedSymbols = [recivedSymbols; x_decoded];
            x_estimate_real = [real(recivedSymbols); imag(recivedSymbols)];
            
            if phase == 2
                estimate_symbols = blackBox(net, x_decoded);
                prediction = [prediction; estimate_symbols];
            end
    end
    
    if phase == 1
        net = fitcnet(x_estimate_real,modData_real,"LayerSizes",[10 10 10])
        figure;
        plot(net.TrainingHistory.TrainingLoss)
    end  
end

%% convert back to bits

        % figure out which symbol maps to which bits.
        outputVector =  convertToBits(modData, data, prediction);

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