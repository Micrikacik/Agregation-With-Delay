function postfix = FPFilePostfix(d, stepDelay)

arguments
    d 
    stepDelay 
end

postfix = sprintf("%iD_%i", d, stepDelay);
