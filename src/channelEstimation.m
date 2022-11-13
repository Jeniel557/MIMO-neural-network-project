function [H_estimate] = channelEstimation(R,t1,t2)

	R_squarMatrix = [R(1) , R(3) ; R(2) , R(4)];
    H_estimate = abs( R_squarMatrix*inv([t1, t2]));
  
end

