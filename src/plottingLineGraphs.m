close all;
clear all;
clc;
format long
% This is the main file, dont run BBS or BBSoneNN


i = sqrt(-1);

trainingTimeNN1 = 0;
trainingTimeNN2 = 0;
testingTime = [];

snrData = zeros(8,2);
% snrData(3,:) = 4
rx = 2;
tx = 2;
H_store = [];
R_store =[];
BER_prediction1 = [];
BER_prediction2 = [];


numBits = 200000;
numSymbols = 4;

symbolSize = log2(numSymbols);

data_train = randi([0 1],numBits,1); % Generate training data.
data = data_train;

figure

SNR = [1:2:20];

for snr = SNR
    
    BER_prediction1 = [];
BER_prediction2 = [];
    
    disp("-----------------------------------------------");
    snr
    rng default
    hAWGN = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (SNR)',...
        'SNR',snr,...
        'SignalPower',1);
    
    modData = [];
    Rayscale = 0.75;
    recivedSymbols = [];
    prediction = [];
    
    H_not_pilot = [];
    TransmissionCounter = 1; % starts as one to account for the first transmission.
    
    xSize = (numBits/symbolSize);
    
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
    X_total_real = [real(X_total); imag(X_total)];
    
    
    
    % add nosie
    T_total = step(hAWGN,T_total);
    
    
    
    
    
    
    T_counter = 1;
    H_counter = 1;
    X_counter = 1;
    
    % disp("-----------------------------------------------")
    
    H_store = zeros((xSize*2),1);
    Hstore_counter = 1;
    R_store =  zeros((xSize*2),1);
    rStore_counter = 1;
    T_bar = zeros((xSize*2),1);
    tBar_counter = 1;
    
    recivedSymbols = zeros((xSize),1);
    recivedSymbols_counter = 1;
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
        x_decoded_hat= decoder(H_actual, R);
        
        % Precoder
        t1_bar = [x1; x2]; % Assume accurate decoder
        t2_bar = [conj(-x2); conj(x1)]; % Assume accurate decoder
        T = [t1_bar; t2_bar];
        
        % Training information storeage
        %     H_store = [H_store; [H_actual(:,1);H_actual(:,2)]]
        H_store(Hstore_counter) =  H_actual(1,1);
        H_store(Hstore_counter+1) =  H_actual(2,1);
        H_store(Hstore_counter+ 2) =  H_actual(1,2);
        H_store(Hstore_counter + 3) =  H_actual(2,2);
        Hstore_counter = Hstore_counter + 4;
        
        %     R_store = [R_store; R]
        R_store(rStore_counter : rStore_counter+3) = R;
        rStore_counter = rStore_counter + 4;
        
        %      recivedSymbols = [recivedSymbols; x_decoded_hat]
        recivedSymbols(recivedSymbols_counter : recivedSymbols_counter + 1) = x_decoded_hat;
        recivedSymbols_counter = recivedSymbols_counter + 2;
        
        T_bar(tBar_counter : tBar_counter + 3) =  T;
        tBar_counter = tBar_counter + 4;
        
    end
    x_estimate_real = [real(recivedSymbols); imag(recivedSymbols)];
    
    
    
%     display("Training nn1 ")
    tNN1 = tic;
    net = fitcnet(x_estimate_real, X_total_real,"LayerSizes",[10 10 10 6]);
    trainingTimeNN1 = toc(tNN1);
%     display(trainingTimeNN1)
    
    input = [real(T_bar) ,imag(T_bar) , real(R_store), imag(R_store)];
%     display("Training nn2 ")
    tnn2 = tic;
    CEnet = fitrnet(input, H_store,"LayerSizes",[10 10 10]);
    trainingTimeNN2 = toc(tnn2);
%     display(trainingTimeNN2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Testing
    
    
    numBits = 80000;
    symbolSize;
    numSymbols;
    xTestSize = (numBits/symbolSize);
    

    traningSNRs = [1:20];
    
    for snrtraining = traningSNRs 
        rng default
        hAWGN = comm.AWGNChannel(...
            'NoiseMethod','Signal to noise ratio (SNR)',...
            'SNR',snrtraining,...
            'SignalPower',1);
        
        
    
        
        data_test = randi([0 1],numBits,1); % generate testing data
        data = data_test;
        
        T_total_test = zeros((xTestSize*2),1);
        t_total_test_counter = 1;
        H_total_test = [];
        X_total_test = zeros((xTestSize),1);
        x_tot_test_counter = 1;
        T1_total_test = zeros((xTestSize),1);
        T2_total_test = zeros((xTestSize),1);
        t1and2_counter = 1;
        modData = [];
        prediction1 = zeros((xTestSize),1);
        pred1and2_counter = 1;
        prediction2 = zeros((xTestSize),1);
        pilotRecalculation = 100;
        
        TransmissionCounter = 0;
        for j = 1:2*symbolSize:length(data)
            
            TransmissionCounter = TransmissionCounter+1;
            
            % Modulate the data into symbols.
            [x1, x2] = modulation(data, j, numSymbols);
            %         X_total_test = [X_total_test; x1 ; x2];
            X_total_test(x_tot_test_counter) = x1;
            X_total_test(x_tot_test_counter + 1) = x2;
            x_tot_test_counter = x_tot_test_counter + 2;
            % Encode the symbols with the Alamoteee
            t1 = [x1; x2];
            t2 = [conj(-x2); conj(x1)];
            T1_total_test(t1and2_counter : t1and2_counter + 1) = t1;
            T2_total_test(t1and2_counter : t1and2_counter + 1) = t2;
            t1and2_counter = t1and2_counter +2;
            
            %         T1_total_test = [T1_total_test; t1];
            %         T2_total_test = [T2_total_test; t2];
            
            % channel matrix
            
            H_actual = KnownChannelMatrixCreation(Rayscale,rx,tx);

            R = [H_actual*t1; (H_actual*t2)];
 
            T_total_test(t_total_test_counter :t_total_test_counter+  3) = R;
            t_total_test_counter = t_total_test_counter + 4;
            
        end
        
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
            prediction1(pred1and2_counter : pred1and2_counter+1) = x_bar;
            prediction2(pred1and2_counter : pred1and2_counter+1) = x_bar_prime;
            
            pred1and2_counter = pred1and2_counter + 2;
            
            
        end
        timerTemp = toc(tTesting);
        testingTime = [testingTime ; [snr , timerTemp ]];
        
        %% Convert back to bits
        
        
        outputVector =  convertToBits(prediction1, numSymbols);
        numberOfBitErrors = sum(abs(outputVector -  data));
        BER = numberOfBitErrors/size(data,1);
        Percent_Accuracy = (1 -BER )*100;
        BER_prediction1 = [BER_prediction1, BER];
        
        outputVector =  convertToBits( prediction2, numSymbols);
        numberOfBitErrors = sum(abs(outputVector -  data));
        BER = numberOfBitErrors/size(data,1);
        Percent_Accuracy = (1 -BER )*100;
        BER_prediction2 = [BER_prediction2, BER];
        
        
    end
    
%     traningSNRs
%     BER_prediction2
    tableVales = [traningSNRs' ,BER_prediction1' , BER_prediction2' ]
    
    semilogy(traningSNRs,BER_prediction2)
    hold on;
end


%legend("1 dB", "2 dB","3 dB"," 4 dB", "5 dB", "6 dB", "7 dB","8 dB"," 9 dB", "10 dB", "11 dB", "12 dB","13 dB"," 14 dB", "15 dB", "16 dB", "17 dB","18 dB","19 dB", "20 dB")
legend("1 dB","3 dB", "5 dB", "7 dB"," 9 dB", "11 dB","13 dB","15 dB", "17 dB","19 dB")

title("BER comparing the iterative and the black box approach");
xlabel("SNR (dB)");
ylabel("Bit Error Rate (BER)");



% tableVales = [SNR' ,BER_prediction1' , BER_prediction2' ]
% display("Time taken to train nn1 = " + trainingTimeNN1);
% display("Time taken to train nn2 = " + trainingTimeNN2);
% % display("Time taken to test " + numBits + " bits is  = " + timerTemp);
%
% sumTime = sum(testingTime(:,2));
% averageTestingTime = sumTime/30
% display("Above is the time taken ot test "+ numBits + " bits ");

