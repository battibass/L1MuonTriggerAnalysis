#!/usr/bin/python

from subprocess import call, Popen, PIPE
import sys

eosCmd = '/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select'
eosBaseDir = "/eos/cms/store/user/battilan/L1Trigger/L1Muon2015/" 


def eosFindDirs( basePath ) :

    eosFindCmd = Popen([eosCmd + ' find -d ' + basePath], shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = eosFindCmd.communicate()

    dirs = []
    
    if eosFindCmd.returncode == 0 :
        for line in stdout.split("\n") :
            if len(line) > 0 : # CB minor hack I want just v4
                dirs.append(line)
    else :
        print "Error in eosFindDirs"

    return dirs


def eosGetFiles ( path ) :

    eosLsCmd = Popen([eosCmd + ' ls ' + path], shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = eosLsCmd.communicate()

    files = []
    
    if eosLsCmd.returncode == 0 :
        for line in stdout.split("\n") :
            if len(line) > 0 :
                files.append(line)
    else :
        print "Error in eosGetFiles"

    return files


def eosCopy ( inputFile, outputPath ) :

    #print outputPath
    #print '/afs/cern.ch/project/eos/installation/0.3.4/bin/eos.select cp ' + inputFile + " " + outputPath
    call([eosCmd + ' cp ' + inputFile + " " + outputPath], shell=True)


def mkDir( path ) :

    call(["mkdir", "-p",  path])

def mergeFiles( outputBaseDir, mergedFileTag, filesToMerge) :

    fileName = outputBaseDir + "/L1Tree_" + mergedFileTag + ".root"

    #print ["hadd", "-f",  fileName] + filesToMerge
    call(["hadd", "-f",  fileName] + filesToMerge)


def copyAndMerge ( inputBasePath, outputBasePath, nFilesToCopy, fileNameTag ) :

    for dir,files in getFilesToCopy( inputBasePath, nFilesToCopy, fileNameTag ).items() :

        outputPath = (outputBasePath + "/" + dir.replace( eosBaseDir, "" )).replace("//","/")

        mergedFileTag = (dir.replace( eosBaseDir, "" )).replace("/","_")
        if mergedFileTag[0] == "_" :
            mergedFileTag = mergedFileTag[1:]
        if mergedFileTag[-1] == "_" :
            mergedFileTag = mergedFileTag[:-1]

        mkDir( outputPath )

        filesToMerge = []
        
        for file in files:
            inputPath = (dir + file).replace("//","/")
            eosCopy( inputPath, outputPath)
            filesToMerge.append(outputPath + "/" + file)

        mergeFiles(outputBaseDir, mergedFileTag, filesToMerge)
        
        

def getFilesToCopy( basePath, nFilesToCopy, fileNameTag ) :

    dirWithFiles = []
    for dir in eosFindDirs( basePath ) :
        nFiles = int(dir.split()[2].strip("nfiles="))
        if nFiles > 0 :
            dirWithFiles.append(dir.split()[0])

    filesToCopy = {}

    for dir in dirWithFiles :
        files = []
        for file in eosGetFiles( dir ) :
            if file.find(fileNameTag) != -1 :
                if nFilesToCopy > len(files) :
                    files.append( file )

        filesToCopy[dir] = files

    return filesToCopy


nArgs = len(sys.argv)
if nArgs < 3 :
    print "Usage: " , sys.argv[0], " input_base_dir(eos) output_base_dir(local) [#_of_files_to_copy(per sample)]"
    sys.exit(1)

inputBaseDir  = sys.argv[1]
outputBaseDir = sys.argv[2]

if nArgs == 4 :
    nFilesToCopy = int(sys.argv[3])
else :
    nFilesToCopy = 99999

print sys.argv

copyAndMerge( inputBaseDir, outputBaseDir, nFilesToCopy, 'L1Tree' )
    
        
    
