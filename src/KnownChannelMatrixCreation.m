function [H] = KnownChannelMatrixCreation(scale,h_rows,h_columns )
%KNOWNCHANNELMATRIXCREATION Summary of this function goes here
%   Detailed explanation goes here
% not filtered, 
h11 = raylrnd(scale);
h22 = raylrnd(scale);
h12 = raylrnd(scale);
h21 = raylrnd(scale);
H = [ h11, h12; h21, h22];


end

