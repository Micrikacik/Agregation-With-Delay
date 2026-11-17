function folderPath = FPFolderPath(d, stepDelay, make)

arguments
    d
    stepDelay
    make = false
end     

folderPath = sprintf("Fokker-Planck/%iD/stepDelay%i", d, stepDelay);
if make && ~exist(folderPath,"dir")
    mkdir(folderPath)
end