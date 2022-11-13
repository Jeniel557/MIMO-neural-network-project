load('tTotalVariable.mat')

numBits = 100000;
numSymbols = 16;
symbolSize = log2(numSymbols);
Rayscale = 0.75;


BER_prediction1 = [];
BER_prediction2 = [];
numBits = 100000;
SNR = [1:2:20];

for snr = SNR
    disp("-----------------------------------------------");
    snr
    rng default
    hAWGN = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR',snr,...
        'SignalPower',1);
    
    
    data_test = randi([0 1],numBits,1); % generate testing data
    data = data_test;
    
    T_total_test = [];
    H_total_test = [];
    X_total_test = [];
    T1_total_test = [];
    T2_total_test = [];
    modData = [];
    prediction1 = [];
    prediction2 = [];
    pilotRecalculation = 100;
    
    TransmissionCounter = 0;
    for j = 1:2*symbolSize:length(data)
        
        TransmissionCounter = TransmissionCounter+1;
        
        % Modulate the data into symbols.
        [x1, x2] = modulation(data, j, numSymbols);
        X_total_test = [X_total_test; x1 ; x2];
        
        % Encode the symbols with the Alamoteee
        t1 = [x1; x2];
        t2 = [conj(-x2); conj(x1)];
        T1_total_test = [T1_total_test; t1];
        T2_total_test = [T2_total_test; t2];
        
        % channel matrix
        
        H_actual = KnownChannelMatrixCreation(Rayscale,rx,tx);
        H_total_test = [H_total_test ; H_actual];
        
        % Recieved values
        R = [H_actual*t1; (H_actual*t2)];
        T_total_test = [T_total_test ; R];
    end
    
    % disp("-----------------------------------------------");
    
    T_total_test = step(hAWGN,T_total_test);
    
    T_counter = 1;
    H_counter = 1;
    X_counter = 1;
    T1_counter = 1;
    T2_counter = 1;
    TransmissionCounter = 0;
    X_decoded_total = [];
    recivedSymbols_test = [];
    
    tTesting = tic;
    for j = 1:2*symbolSize:length(data)
        TransmissionCounter = TransmissionCounter+1;
        
        % reconstruct R and H_est_ls and X
        R = [T_total_test(T_counter); T_total_test(T_counter + 1);T_total_test(T_counter +2);T_total_test(T_counter+3)];
        T_counter = T_counter + 4;
        
        t1 = [T1_total_test(T1_counter) ;T1_total_test(T1_counter+1) ];
        t2 = [T2_total_test(T1_counter) ;T2_total_test(T1_counter+1) ];
        T1_counter = T1_counter+ 2;
        
        % recalculate the pilot and save it
        if TransmissionCounter == 1
            H_est_ls = channelEstimation(R,t1,t2);
        end
        % set the H to be recalculated in the next transmision
        if TransmissionCounter == pilotRecalculation
            TransmissionCounter = 0;
        end
        
        % Decoder
        R(3:4) = conj(R(3:4));
        x_decoded_hat= decoder(H_actual, R);
        
        % First network
        pred = predict(net,[real(x_decoded_hat); imag(x_decoded_hat)]);
        x_bar = [pred(1) + pred(3)*i; pred(2) + pred(4)*i];
        
        %Back into the precoder to produce T.
        t1_bar = [x_bar(1); x_bar(2)]; % Assume accurate decoder
        t2_bar = [conj(-x_bar(2)); conj(x_bar(1))]; % Assume accurate decoder
        T = [t1_bar; t2_bar];
        
        % Second network to output H_bar
        input = [real(T) ,imag(T) , real(R), imag(R)];
        H_bar = predict(CEnet, input); % Assume accurate prediction
        H_bar = [H_bar(1:2), H_bar(3:4)];
        
        % Decoder from second network output.
        x_decoded_bar_prime = decoder(H_bar, R);
        
        % First network Again, final prediction
        pred = predict(net,[real(x_decoded_bar_prime); imag(x_decoded_bar_prime)]);
        x_bar_prime = [pred(1) + pred(3)*i; pred(2) + pred(4)*i];
        
        % Training information storeage
        prediction1 = [prediction1; x_bar]; % storage
        prediction2 = [prediction2; x_bar_prime];
        
        
        modData = [modData; [x1; x2]];
        modData_real = [real(modData); imag(modData)];
        
        recivedSymbols_test = [recivedSymbols_test; x_decoded_hat];
        X_decoded_total = [real(recivedSymbols_test); imag(recivedSymbols_test)];
        
    end
    timerTemp = toc(tTesting)
    testingTime = [testingTime ; [snr , timerTemp ]];
    
    %% Convert back to bits
    
    
    
    outputVector =  convertToBits(prediction1, numSymbols);
    numberOfBitErrors = sum(abs(outputVector -  data));
    BER = numberOfBitErrors/size(data,1);
    Percent_Accuracy = (1 -BER )*100
    BER_prediction1 = [BER_prediction1, BER];
    
    outputVector =  convertToBits( prediction2, numSymbols);
    numberOfBitErrors = sum(abs(outputVector -  data));
    BER = numberOfBitErrors/size(data,1);
    Percent_Accuracy = (1 -BER )*100
    BER_prediction2 = [BER_prediction2, BER];
end







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
plot(prediction2,'ob')
xlim([-2 2])
ylim([-2 2])
xlabel("Real component")
ylabel("Imaginary component")
legend("Transmitted", "recieved")
title("Decoded vs transmitted symbol mapping")


figure
semilogy(SNR,BER_prediction1)
hold on;
semilogy(SNR,BER_prediction2)
legend("Black Box", "Iterative")
title("BER comparing the iterative and the black box approach");
xlabel("SNR (dB)");
ylabel("Bit Error Rate (BER)");






tableVales = [SNR' ,BER_prediction1' , BER_prediction2' ]
display("Time taken to train nn1 = " + trainingTimeNN1);
display("Time taken to train nn2 = " + trainingTimeNN2);
% display("Time taken to test " + numBits + " bits is  = " + timerTemp);
testingTime
sumTime = sum(testingTime(:,2));
averageTestingTime = sumTime/30
display("Above is the time taken ot test "+ numBits + " bits ");