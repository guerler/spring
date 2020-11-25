def createFile(identifier, databaseIndex, database, outputName):
    start = -1
    size = 0
    with open(databaseIndex) as file:
        for line in file:
            cols = line.split()
            if identifier == cols[0]:
                start = int(cols[1])
                size = int(cols[2])
                break
    if start != -1 and size > 0:
        with open(database) as file:
            file.seek(start)
            content = file.read(size)
            outputFile = open(outputName, "w")
            outputFile.write(content)
            outputFile.close()
        return True
    else:
        return False
