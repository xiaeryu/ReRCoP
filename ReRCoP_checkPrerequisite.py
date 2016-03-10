import os
import re
import imp
import commands


def checkModule(module):
    '''
    This function raises an error if a python module is not installed.
    
    @param module: the module name
    '''
    try:
        imp.find_module(module)
    except:
        raise ImportError("Please check the installation of the python module '%s'!" % module)


def checkCommand(command):
    '''
    This function raises an error if a command cannot be executed.

    @param command: the command name
    '''
    status, result = commands.getstatusoutput(command)
    if status != 0:
        raise OSError("Cannot execute the command: '%s'!" % command)


def checkFile(fileName):
    '''
    This function raises an error if a file does not exist.

    @param fileName: the file name to check
    '''
    if not os.path.isfile(fileName):
        raise IOError("Invalid input file: '%s'!" % fileName)


def checkDir(dirName):
    '''
    This function creates a directory if there is not one and raises an error if
    cannot create the directory.

    @param dirName: the directory name
    '''
    if not os.path.isdir(dirName):
        try:
            os.makedirs(dirName)
        except:
            raise OSError("Cannot create directory: '%s'!" % dirName)


def checkName(inFasta):
    '''
    This function checks the name of sequences in a fasta file and raises
    an error if there are duplicated names.

    @param inFasta: input fasta file
    '''
    storage = set()
    duplicate = set()
    pat = re.compile(">(\S+)")
    inH = open(inFasta)
    for line in inH:
        if line.startswith('>'):
            name = re.search(pat,line).group(1)
            if name in storage and name not in duplicate:
                duplicate.add(name)
            else:
                storage.add(name)
    inH.close()
    if len(duplicate) > 0:
        raise IOError("Please check the sequence names in '%s'. Duplicate names: %s" % (inFasta, str(duplicate)))


def checkLen(inFasta):
    '''
    This function checks the lengths of the fasta sequences and raises
    an error if the lengths are not equal.

    @param: inFasta: a fasta object
    '''
    length = len(inFasta[inFasta.keys()[0]])
    for seq in inFasta:
        if len(inFasta[seq]) != length:
            raise IOError("Different fasta sequence lengths! Please check!")
