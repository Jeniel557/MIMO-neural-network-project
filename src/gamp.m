function [x] = gamp(y,H, noise_var, numIterations)

    u_x = 1;
    x = 0;
    s = 0;
    for t = 1:numIterations
        for i = 1:length(y)
            
            u_p = sum((H(i,:).^2)*u_x);

            p = sum( H(i,:) * x - u_p * s);

            s = (y - p)/(u_p + noise_var);

            u_s = 1/(u_p + noise_var);

            u_r = inv(sum((H(i,:).^2)*u_s));

            r = x + u_r * sum( H(i,:) * s);

            x = G(r, 0.6);
            
            if i < 5
               u_x = u_r * sum(r)/length(r);
            else
                u_x = u_r * sum(r(i - 5: i))/5;
            end
        end
    end
end

