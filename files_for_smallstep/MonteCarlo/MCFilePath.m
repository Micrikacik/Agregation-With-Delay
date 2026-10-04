function filePath = MCFilePath(folderPath, fileName, postfix)

filePath = sprintf("%s/%s_%s.mat", folderPath, fileName, postfix);