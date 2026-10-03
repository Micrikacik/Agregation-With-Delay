function postfix = MCGroupFilePostfix(group)

arguments
    group (1,1) double {mustBeInteger, mustBePositive}
end

postfix = sprintf("group%i", group);