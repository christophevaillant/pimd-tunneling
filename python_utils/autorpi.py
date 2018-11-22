from tempfile import mkstemp
import os
import numpy as np
import shutil
from subprocess import call

def replace(file_path, pattern, subst):
    "Replaces a string pattern with subst in file file_path."
    #Create temp file
    fh, abs_path = mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    os.close(fh)
    #Remove original file
    os.remove(file_path)
    #Move new file
    shutil.move(abs_path, file_path)


#---------------------------------------------------------
#Work out spacings and generate upper and lower bounds
lowestbeta= 66000.0#lowest beta
highestbeta= 75000.0 #corresponding highest beta

factr= 1.1

numbeta= int(round(np.log(highestbeta/lowestbeta)/np.log(factr)))+1

print numn, numbeta

betaranges=np.logspace(0.0, numbeta, endpoint=False,num=numbeta,base=factr)*lowestbeta

nranges=np.logspace(0.0, numn, endpoint=False,num=numn, base=factr)*lowestn
print betaranges
print nranges
#---------------------------------------------------------
#Loop over resonances and submit molscat jobs to locate resonances.
#Assumes the input script is properly set up and requires only changing
#field ranges

inputfilename= "rpinml.in" #filename of molscat input script to be copied
scriptfilename= "rpi_A.sh" #filename of submissions script

#Set up the directory structure:
dircontents= os.listdir(".")
workingdir="BetaN"
if workingdir in dircontents:
    os.chdir(workingdir)
    newdircontents= os.listdir(".")
    for direct in newdircontents:
        os.chdir(direct)
        filelist= os.listdir(".")
        for fil in filelist:
            os.remove(fil)
        os.chdir("..")
        os.rmdir(direct)
else:
    os.mkdir(workingdir)
    os.chdir(workingdir)


#Strings to be matched in the molscat input file for field range:
#NOTE THAT THESE STRINGS NEED TO MATCH THE INPUT FILE EXACTLY!
nstring= "  n="
betastring= "  beta="

#String to go in submission script to actually run molscat:
scriptstring= "../../../../../build/rpi_watdim_ser < "

#template directory name:
templatename= "beta" 
memorystring= "##  #SBATCH --mem=20480"

#Recursively prepare directories and submit scripts
for i in range(len(betaranges)):
    for j in range(len(nranges)):
        # print betaranges[i], int(nranges[j])
        dirname=templatename+str(i) + "_n" + str(j)
        infilename= templatename+str(i)+"_n"+str(j) + ".in"
        outfilename= templatename+str(i) +"_n"+str(j) + ".out"
        os.mkdir(dirname)
        os.chdir(dirname)
        os.symlink("../../well1.dat","./well1.dat")
        os.symlink("../../well2.dat","./well2.dat")
        os.symlink("../../masses.dat","./masses.dat")
        os.symlink("../../path.xyz","./path.xyz")
        shutil.copyfile("../../"+ inputfilename, infilename)
        shutil.copyfile("../../"+ scriptfilename, scriptfilename)
        replace(infilename,nstring, nstring+ str(int(nranges[j])) + ",")
        replace(infilename,betastring, betastring+ str(betaranges[i]) + "D0,")
        replace(scriptfilename,scriptstring, scriptstring + infilename + "> " + outfilename)
        # if (nranges[j] > 700):
        #     replace(scriptfilename,memorystring, memorystring[4:])
        call("sbatch " + scriptfilename,shell=True)
        os.chdir("..")
