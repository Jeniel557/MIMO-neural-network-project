close all;
clear all;
clc;
% This is the main file, dont run BBS or BBSoneNN


i = sqrt(-1);

numDiffererntSNRs = 6;
snrData = zeros(8,2);
% snrData(3,:) = 4
rx = 2;
tx = 2;

numBits = 100;
numSymbols = 4;
symbolSize = log2(numSymbols);

data_train = randi([0 1],numBits,1); % Generate training data.
data = data_train;

recivedSymbols = [];
estimatedSymbols =[];
H_store = [];
R_store =[];
T_bar = [];
TransmissionCounter = 1; % starts as one to account for the first transmission.
T_total = zeros(numBits/symbolSize,1);
H_total = zeros(numBits/symbolSize,2);
X_total = zeros(numBits/symbolSize,1);
X_estimate = zeros(numBits/symbolSize,1);

SNR = (1:40);
BER_store = [];
for snr = SNR
    disp("------------- " + snr +" --------------")
    counter = 1;
    for j = 1:2*symbolSize:length(data)

        % Modulate the data into symbols.
        [x1, x2] = modulation(data, j, numSymbols);

        % Encode the symbols with the Alamoteee
        t1 = [x1; x2];

        % channel matrix
        scale = 0.5;
        H_actual = KnownChannelMatrixCreation(scale,rx,tx);
        H_total(counter, :) = H_actual(1,:);

        % Recieved values
        R_no_noise = [H_actual*t1];

        T_total((counter:counter+1), 1) = R_no_noise;
        X_total((counter:counter+1), 1) = [x1; x2];
        H_total(counter + 1, :) = H_actual(2,:);

        counter = counter + 2;
    end

    rng default
    hAWGN = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR',snr,...
        'SignalPower',1);

    X_total_real = [real(X_total); imag(X_total)];
    T_total = step(hAWGN,T_total);

    counter2 = 1;
    X_decoded = zeros(numBits/symbolSize,1);

    for z = 1:2*symbolSize:length(data)
        % construct R and H_est_ls and X
        R = [T_total(counter2); T_total(counter2+ 1)];

        H_actual = [H_total(counter2,:); H_total(counter2+1, :)];
        
        % Decoder
        x_d = inv(H_actual)*R;

        % Argmax
        x_estimated = [0;0];
        x_estimated(1) = argmax(x_d(1), numSymbols);
        x_estimated(2) = argmax(x_d(2), numSymbols);

        X_decoded((counter2: counter2 + 1),1) = x_d;
        X_estimate((counter2: counter2 + 1),1) = x_estimated;
        counter2 = counter2 + 2;
    end
     
       outputVector = convertToBits(X_estimate, numSymbols);
    numberOfBitErrors = sum(abs(outputVector -  data));
    BER = numberOfBitErrors/size(data,1)
    Percent_Accuracy = (1 -BER )*100;
        
    BER_store = [BER_store;BER];
 
end




figure;
plot(X_total, 'xr','MarkerSize',15)
hold on;
plot(X_decoded,'ob')
xlim([-2 2])
ylim([-2 2])
xlabel("Real component")
ylabel("Imaginary component")
legend("Transmitted", "recieved")
%title("Decoded vs transmitted symbol mapping")


figure;
plot(X_total, 'xr','MarkerSize',15)
hold on;
plot(X_estimate,'ob')
xlim([-2 2])
ylim([-2 2])
xlabel("Real component")
ylabel("Imaginary component")
legend("Transmitted", "recieved")
%title("Decoded vs transmitted symbol mapping")

figure
semilogy(SNR,BER_store)
legend("Black Box")
title("BER comparing the iterative and the black box approach");
xlabel("SNR (dB)");
ylabel("Bit Error Rate (BER)");





