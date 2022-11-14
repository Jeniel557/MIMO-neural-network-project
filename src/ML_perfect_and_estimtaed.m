close all;
clear all;
clc;
% This is the main file for Maximun likeyhood with perfect H. Run this
% first.


i = sqrt(-1);

numDiffererntSNRs = 6;
snrData = zeros(8,2);

rx = 2;
tx = 2;

numBits = 1000000;
numSymbols = 4;
Rayscale = 0.5;

symbolSize = log2(numSymbols);
xSize = (numBits/symbolSize);

SNR = [1:2:40];
tableData = zeros(length(SNR), 3);
tableDataCounter = 1;
for snr = SNR
    snr
 
    rng default
    hAWGN = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR',snr,...
        'SignalPower',1);
    data_train = randi([0 1],numBits,1); % Generate training data.
    data = data_train;

    recivedSymbols = [];
    estimatedSymbols_estiamted =zeros( xSize,1);
    estimatedSymbols_perfect= zeros( xSize,1);
    estimated_symbol_counter = 1;
  
%     R_store =[];
  
T_total = zeros((xSize*2),1);
H_total = zeros( xSize,2);
X_total = zeros( xSize,1);
t_counter = 1;
x_counter = 1;
h_counter = 1;



for j = 1:2*symbolSize:length(data)
    
    % Modulate the data into symbols.
    [x1, x2] = modulation(data, j, numSymbols);
    X_total(x_counter) = x1;
    X_total(x_counter+1) = x2;
    x_counter = x_counter + 2;
    
    % Encode the symbols with the Alamoteee
    t1 = [x1; x2];
    t2 = [conj(-x2); conj(x1)];
    
    % channel matrix
    
    H_actual = KnownChannelMatrixCreation(Rayscale,rx,tx);
    H_total( h_counter: h_counter + 1, 1:2) = H_actual;
    h_counter  = h_counter + 2;
    
    % Recieved values
    R_no_noise = [H_actual*t1; (H_actual*t2)];
    T_total(t_counter) = R_no_noise(1);
    T_total(t_counter + 1) = R_no_noise(2);
    T_total(t_counter + 2) = R_no_noise(3);
    T_total(t_counter + 3) = R_no_noise(4);
    t_counter = t_counter + 4;
end
    %% adding noise
    T_total = step(hAWGN,T_total);

    T_counter = 1;
    H_counter = 1;
    X_counter = 1;

    %% argmax decoding
    TransmissionCounter = 0;
    pilotRecalculation = 1;
    numTimesRecalc = 0;
    for j = 1:2*symbolSize:length(data)
        TransmissionCounter = TransmissionCounter+1;


        % construct R and H_est_ls and X
        R = [T_total(T_counter); T_total(T_counter + 1);T_total(T_counter +2);T_total(T_counter+3)];
        T_counter = T_counter + 4;

        H_actual = [H_total(H_counter, :)  ;  H_total(H_counter+1, :)  ];
        H_counter = H_counter+2;

        % recalculate the pilot and save it
        if TransmissionCounter == 1
            H_est_ls = channelEstimation(R,t1,t2);
            numTimesRecalc = numTimesRecalc + 1;
        end
        % set the H to be recalculated in the next transmision
        if TransmissionCounter == pilotRecalculation
            TransmissionCounter = 0;
        end


        x1 = X_total(X_counter);
        x2 = X_total(X_counter+1);
        X_counter = X_counter+2;

        % Decoder
        R(3:4) = conj(R(3:4));
        x_decoded_estimate= decoder(H_est_ls, R);

        % Argmax
        x_estimated = [0;0];
        x_estimated(1) = argmax(x_decoded_estimate(1));
        x_estimated(2) = argmax(x_decoded_estimate(2));

        % Training information storeage
        estimatedSymbols_estiamted(estimated_symbol_counter : estimated_symbol_counter +1) = x_estimated;

        % perfect decoding
        x_decoded_perfect= decoder(H_actual, R);
        % Argmax
        x_estimated = [0;0];
        x_estimated(1) = argmax(x_decoded_perfect(1));
        x_estimated(2) = argmax(x_decoded_perfect(2));


        estimatedSymbols_perfect(estimated_symbol_counter : estimated_symbol_counter +1) = x_estimated;
        estimated_symbol_counter = estimated_symbol_counter + 2;
    end


    %% Convert back to bits
    outputVector = convertToBits(estimatedSymbols_estiamted);
    numberOfBitErrors = sum(abs(outputVector -  data));
    BER_estimatted = numberOfBitErrors/size(data,1)
    Percent_Accuracy = (1 -BER_estimatted )*100;


    outputVector = convertToBits(estimatedSymbols_perfect);
    numberOfBitErrors = sum(abs(outputVector -  data));
    BER_perfect = numberOfBitErrors/size(data,1)
    Percent_Accuracy = (1 -BER_perfect )*100;

    tableData(tableDataCounter,:) = [snr; BER_estimatted; BER_perfect];
    tableDataCounter = tableDataCounter + 1;

end

tableData
display("Total times recalcualted is " + numTimesRecalc);

figure;
plot(X_total, 'xr','MarkerSize',15)
hold on;
plot(recivedSymbols,'ob')
xlim([-2 2])
ylim([-2 2])
xlabel("Real component")
ylabel("Imaginary component")
legend("Transmitted", "recieved")
title("Decoded vs transmitted symbol mapping")

figure;
plot(X_total, 'xr','MarkerSize',15)
hold on;
plot(estimatedSymbols_estiamted,'ob')
xlim([-2 2])
ylim([-2 2])
xlabel("Real component")
ylabel("Imaginary component")
legend("Transmitted", "recieved")
title("Decoded vs transmitted symbol mapping")


figure
semilogy(tableData(:,1), tableData(:,2))
hold on
semilogy(tableData(:,1), tableData(:,3))
legend("ML with estimated H", "ML with perfect H ")
title("BER comparing the iterative and the black box approach");
xlabel("SNR (dB)");
ylabel("Bit Error Rate (BER)");

