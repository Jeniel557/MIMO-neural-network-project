function [H_estimate] = channelEstimation(R,t1,t2)
%CHANNELESTIMATION Summary of this function goes here
%   Detailed explanation goes here

R_squarMatrix = [R(1) , R(3) ; R(2) , R(4)];
        symbolsToTransmis = [t1 , t2];
%         display(H)
        H_estimate = abs( R_squarMatrix*inv(symbolsToTransmis));


end

