function folderPath = MCFolderPath(delayType, d, stepDelay, make)

arguments
    delayType
    d
    stepDelay
    make = false
end

folderPath0 = sprintf("MonteCarlo");
if make && ~exist(folderPath0,"dir")
    mkdir(folderPath0)
end

folderPath1 = sprintf("%s/delayType%s", folderPath0, delayType);
if make && ~exist(folderPath1,"dir")
    mkdir(folderPath1)
end

folderPath2 = sprintf("%s/%iD", folderPath1, d);
if make && ~exist(folderPath2,"dir")
    mkdir(folderPath2)
end        

folderPath3 = sprintf("%s/stepDelay%i", folderPath2, stepDelay);
if make && ~exist(folderPath3,"dir")
    mkdir(folderPath3)
end

folderPath = folderPath3;