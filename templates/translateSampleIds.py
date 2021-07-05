#
#   Translate a file column with a dictionary file
#   ==============================================
#
#   translates the values in a selected column of an input data file
#   to new values, where the translations are provided by a second 
#   'dictionary' file. For the dictionary file, we assume the file has
#   (at least) two columns: entries in the first are the keys, and 
#   entries in the second column are the values.
#
#   General use:
#   $ python translate.py <data file> <target column> <dictionary file>
#
#   The separator/deliminator need not be a space, but it does need to 
#   be the same for both the data file and the dictionary file (and will
#   be the same for the output). The output is terminal standard out, 
#   so you might want to redirect it to a file, i.e. with '>'.
#
#   To get the dictionary for an Illumina sample report, just use:
#
#   $ cat <illumina sample report csv> | tr ',' ' ' | tr '_' ' ' | awk '{print $4" "$2"_"$3"_"$4}' - > <dictionary tsv>
#
#   then call this python script twice to translate the clinical fam columns:
#   $ python translate.py <clinical fam> 0 <dictionary tsv> > temp.fam
#   $ python translate.py temp.fam 1 <dictionary tsv> > <mapped clinical fam>
#
#######################################################################

import sys

SEPARATOR = " "

def main():
    dataFile, targetColumn, dictionaryFile = getCommandLineArguments()
    dictionary = buildDictionary(dictionaryFile)
    translateColumnInDataFile(
        dataFile,
        targetColumn,
        dictionary)

def getCommandLineArguments():
    dataFile = sys.argv[1]
    targetColumn = int(sys.argv[2])
    dictionaryFile = sys.argv[3]
    return [dataFile, targetColumn, dictionaryFile]


def buildDictionary(inputFilePath):
    inputStream = open(inputFilePath, "r")
    dictObject = {}
    lines = inputStream.readlines()
    for line in lines:
        lineArray = line.strip().split(SEPARATOR)
        dictObject[lineArray[0]] = lineArray[1]
    inputStream.close()
    return dictObject


def translateColumnInDataFile(dataFile, targetColumn, dictionary):
    dataStream = open(dataFile, "r")
    dataLines = dataStream.readlines()
    for line in dataLines:
        lineArray = line.strip().split(" ")
        if lineArray[targetColumn] in dictionary.keys():
            lineArray[targetColumn] = dictionary[lineArray[targetColumn]]
        newLine = SEPARATOR.join(lineArray)
        print(newLine)       

if __name__ == '__main__':
    main()