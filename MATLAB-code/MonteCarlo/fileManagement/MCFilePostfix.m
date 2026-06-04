function postfix = MCFilePostfix(delayType, d, stepDelay, group, groupCount)

arguments
    delayType 
    d 
    stepDelay 
    group = -1
    groupCount = -1
end

if (groupCount ~= -1) && (group ~= -1)
    postfix = sprintf("%s_%iD_%i_group%iof%i", delayType, d, stepDelay, group, groupCount);
else
    postfix = sprintf("%s_%iD_%i", delayType, d, stepDelay);
end