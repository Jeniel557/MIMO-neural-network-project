function [output] = G(input, radius)
    
    if input >= 0
        
        if abs(input - 1/sqrt(2)) >= radius
            
            if input > 1/sqrt(2)
                output = 1/sqrt(2) + radius;
            else
                output = 1/sqrt(2) - radius;
            end
        else
            output = input;
        end
    else
        
        if abs(input + 1/sqrt(2)) >= radius
            
            if input > -1/sqrt(2)
                output = -1/sqrt(2) + radius;
            else
                output = -1/sqrt(2) - radius;
            end
        else
            output = input;
        end
        
    end
end

