# INVESTIGATION OF MACHINE LEARNING TECHNIQUES TO IMPROVE MIMO TRANSMISSIONS USING A 2X2 ALAMOUTI STBC


The use of neural networks to replace the traditional means of MIMO signal detection is investigated. A 2x2 MIMO Alamouti STBC communication channel was simulated on MATLAB. A multi layer perceptron (MLP) neural network classifies the received symbols. A second MLP is used to improve channel estimation using an iterative approach. Results show that both these approaches outperform traditional STBC Alamouti decoding, as well as non-Alamouti MIMO communication systems. The Iterative approach achieves a BER of $10^{-5}$ at 25dB, outperforming previous methods.

The src folder contains the implementation code. 

main.m contains the code to run both the Black Box and Iterative approach. 
Adjust the numBits variable in both the training and testing section to control the number of bits to simulate. 

The numSymbols variable changes the number of symbols to be transmitted. 


ThreeIterations.m contains the code to run the Iterative approach with three iterations. 

ML_perfect_and_estimtaed.m contains the code for tradition MIMO signal detection using Maximum Likelihood detection. 
This file contains two methods. One uses a perfect channel estimation, and the other uses a pilot symbol to estimate 
the channel, recalculated ever 10 transmission, controlled by the pilotRecalculation variable.


The non-Alamouti code can be found in the non-Alamouti folder under src
