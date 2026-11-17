function postfix = MCFilePostfix(delayType, d, stepDelay)

arguments
    delayType 
    d 
    stepDelay
end

postfix = sprintf("%s_%iD_%i", delayType, d, stepDelay);