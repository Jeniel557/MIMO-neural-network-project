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
H_store = [];
R_store =[];

SNR = 6;
rng default
hAWGN = comm.AWGNChannel(...
    'NoiseMethod','Signal to noise ratio (SNR)',...
    'SNR',SNR,...
    'SignalPower',1);


numBits = 8000;
numSymbols = 16;
symbolSize = log2(numSymbols);

data_train = randi([0 1],numBits,1); % Generate training data.
data = data_train;

recivedSymbols = [];
estimatedSymbols =[];
H_store = [];
R_store =[];
T_bar = [];
TransmissionCounter = 1; % starts as one to account for the first transmission.
T_total = [];
H_total = [];
X_total = [];

for j = 1:2*symbolSize:length(data)
    
    % Modulate the data into symbols.
    [x1, x2] = modulation(data, j, numSymbols);
    X_total = [X_total; x1 ; x2];
    
    % Encode the symbols with the Alamoteee
    t1 = [x1; x2];
    t2 = [conj(-x2); conj(x1)];
    
    % channel matrix
    scale = 0.75;
    H_actual = KnownChannelMatrixCreation(scale,rx,tx);
    H_total = [H_total ; H_actual];
    
    % Recieved values
    R_no_noise = [H_actual*t1; (H_actual*t2)];
    T_total = [T_total ; R_no_noise];   
end

X_total_real = [real(X_total); imag(X_total)];
T_total = step(hAWGN,T_total);

T_counter = 1;
H_counter = 1;
X_counter = 1;

disp("-----------------------------------------------")

for j = 1:2*symbolSize:length(data)
    
    % construct R and H_est_ls and X
    R = [T_total(T_counter); T_total(T_counter + 1);T_total(T_counter +2);T_total(T_counter+3)];
    T_counter = T_counter + 4;
    
    H_actual = [H_total(H_counter, :)  ;  H_total(H_counter+1, :)  ];
    H_counter = H_counter+2;
    
    x1 = X_total(X_counter);
    x2 = X_total(X_counter+1);
    X_counter = X_counter+2;
    
    % Decoder
    R(3:4) = conj(R(3:4));
    x_decoded= decoder(H_actual, R);
    
    % Argmax
    x_estimated = [0;0];
    x_estimated(1) = argmax(x_decoded(1), numSymbols);
    x_estimated(2) = argmax(x_decoded(2), numSymbols);

    % Training information storeage
    H_store = [H_store; [H_actual(:,1);H_actual(:,2)]];
    R_store = [R_store; R];
    recivedSymbols = [recivedSymbols; x_decoded];
    estimatedSymbols = [estimatedSymbols; x_estimated];
    x_estimate_real = [real(recivedSymbols); imag(recivedSymbols)];
end
sum(estimatedSymbols - X_total)

%% Convert back to bits
outputVector = convertToBits(estimatedSymbols, numSymbols);
numberOfBitErrors = sum(abs(outputVector -  data))
BER = numberOfBitErrors/size(data,1)
Percent_Accuracy = (1 -BER )*100



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
plot(estimatedSymbols,'ob')
xlim([-2 2])
ylim([-2 2])
xlabel("Real component")
ylabel("Imaginary component")
legend("Transmitted", "recieved")
title("Decoded vs transmitted symbol mapping")