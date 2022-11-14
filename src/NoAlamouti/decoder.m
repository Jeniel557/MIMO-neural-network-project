function [x_decoded] = decoder(H_est, R)

    H_bar = [conj(H_est(1,2)), -conj(H_est(1,1)); conj(H_est(2,2)), -conj(H_est(2,1))];
    decoding_H = [H_est; H_bar];
    x_decoded = (inv(decoding_H' * decoding_H)* decoding_H') * R;
    % x_decoded(2) = conj(x_decoded(2));
end

