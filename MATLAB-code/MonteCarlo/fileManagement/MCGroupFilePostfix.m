function postfix = MCGroupFilePostfix(group, groupCount)

arguments
    group (1,1) double {mustBeInteger, mustBePositive}
    groupCount (1,1) double {mustBeInteger, mustBePositive}
end

postfix = sprintf("group%iof%i", group, groupCount);