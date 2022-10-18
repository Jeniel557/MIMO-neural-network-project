function [estimate_symbols] = blackBox(net,x_decoded)

	x_decoded(2) = conj(x_decoded(2));
    x_estimate_real = [real(x_decoded); imag(x_decoded)];
    estimate = predict(net,x_estimate_real);
    estimate_symbols = [estimate(1) + estimate(3)*i; estimate(2) + estimate(4)*i];
end

