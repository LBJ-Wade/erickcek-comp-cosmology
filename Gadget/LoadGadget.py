
import sys
import os

import SetVariables
import MakeDirectories
import RunICs
import RunGadget



##################################################################################
##################################################################################
#################################Data Class#######################################
##################################################################################
##################################################################################
class GadgetData:

    def __init__(self):
        self.variables = {'box_size':'NOT_SET',         #Box size in Mpc/h
                     'n':'NOT_SET',                     #Cube-root of particle number
                     'seed':'13579',                    #Random seed  
                     'power':'power',                   #Name of power spectrum file (without directory or .dat)
                     'gas_flag':'0',                    #0 = no gas; 1 = adiabatic; 2 = chemistry + cooling
                     'nodes':'NOT_SET',                 #Number of nodes to be used
                     'out_code':'NOT_SET',              #File which contains information on when to take snapshots. Full file name is outputs_ + out_code + .txt 
                     'time_limit':'48:00:00',           #Time limit
                     'interconnect':':',                #Interconnect (qdr by default. other allowed value is ddr) ???
                     'mpi_version':'openmpi',           #MPI version (openmpi by default. other allowed values are "intelmpi" or "mpich2")
                     'z_initial':'199',                 #Initial redshift
                     'z_final':'0',                     #Final redshift
                     'delta_file':'NOT_SET',            #Location of input delta file
                     'growth_ratio':'NOT_SET',          #Ratio of D(z_final)/D(z_initial)
                     'chem_flag':'0',                   #???
                     'number_of_snapshots':'10',        #Number of snapshots   
                     'omega_m':'.3061',                   #omega_m
                     'omega_l':'.6939',                   #omega_l
                     'omega_b':'.048093',                  #omega_b
                     'hubble_param':'.6794',              #hubble parameter
                     'restart':'0'}                     #Restart flag

        self.directories = {'root':'NOT_SET',           #Root directory. Where mainDir will be put, and therefore all output
                       'gadget_src':'NOT_SET',          #Contains gadget executable and dummy parameter file
                       'ic_codes':'NOT_SET',            #Location of delta2p and unigrid executables
		       'gen_file' : 'NOT_SET',          #Location of scale factors used for snapshots
                       'input':'NOT_SET'}               #Location of all user input files (power spectrum, outputs file (if needed), 
                                                        #directories file, and variables file)

        self.linearGrowthSpacing = False                #True if a growth_ratio is provided
        self.makeDirectories = False                    #Should directories be made?
        self.runICs = False                             #Should initial conditions be created?
        self.runGadget = False                          #Should Gadget be run?
        self.progressRefreshRate = 0                    #If Gadget is run, how often the user should be updated on its progress (in seconds).
        self.progressFile = 'NOT_SET'
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################



##################################################################################
##################################################################################
#################################MAIN#############################################
##################################################################################
##################################################################################
def go():
    gd = GadgetData()
    checkInput()
    gd.directories['input'] = sys.argv[1]
    initializeVariables(gd.directories['input'] + '/variables.txt',gd.variables)
    if(gd.variables['growth_ratio'] != 'NOT_SET'):
        gd.linearGrowthSpacing = True
    initializeDirectories(gd.directories['input'] + '/directories.txt',gd.directories)
    setRunMode(gd)
    os.chdir(gd.directories['root'])
    setMainDir(gd)
    
    if(gd.makeDirectories):
        MakeDirectories.go(gd)
    if(gd.runICs):
        RunICs.go(gd)
    if(gd.runGadget):
        RunGadget.go(gd)
        
    print("\n\nFinished\n")
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################



##################################################################################
##################################################################################
################################FUNCTIONS#########################################
##################################################################################
##################################################################################
#Makes sure the user gave the program the correct number of arguments.
def checkInput():
    if(len(sys.argv) < 5):
        print """
        LoadGadget takes 5 arguments:
        
            -Directory containing input files
            -Setup directories? (0 or 1)
            -Generate initial conditions? (0 or 1)
            -Run gadget? (0 or 1)
            -(optional) Progress refresh rate in seconds
            -(optional) File to write progress to
        """
        exit(1)
        

#Gets variables from input file and ensures key values are set.
def initializeVariables(variablesFile,variables):
    SetVariables.setVariables(variablesFile,variables)
    if(variables['box_size'] == 'NOT_SET' or variables['n'] == 'NOT_SET' or variables['nodes'] == 'NOT_SET'):
        print('Missing either boxsize, particle number, or node number')
        print('Values: ' + variables['box_size'] + ', ' + variables['n']
              + ', ' + variables['nodes'])
        exit(1)
        

#Gets directory names from input file and checks that all are specified.
def initializeDirectories(directoriesFile,directories):
    SetVariables.setVariables(directoriesFile,directories)
    for key in directories.keys():
        if(directories[key] == -1):
            print('Missing directory \'' + str(key) + '\'')
            exit(1)
        
 
#Sets three boolean variables based on arguments passed to the program which determine what the program will run.
#Also sets a refresh rate for progress updates if Gadget is run.     
def setRunMode(gd):
    if(sys.argv[2] == '1'):
        gd.makeDirectories = True
    if(sys.argv[3] == '1'):
        gd.runICs = True
    if(sys.argv[4] == '1'):
        gd.runGadget = True
    try:
        gd.progressRefreshRate = float(sys.argv[5])
        gd.progressFile = sys.argv[6]
    except:
        pass
        

#Based on gas_flag chooses the main directory name and creates the directory if it does not already exist.   
def setMainDir(gd):
    if(gd.variables['gas_flag'] ==  '0'):
        gd.mainDir = 'gadg_' + gd.variables['n'] + '_' + gd.variables['box_size'] + '_dm'
    elif(gd.variables['gas_flag'] == '1'):
        gd.mainDir = 'gadg_' + gd.variables['n'] + '_' + gd.variables['box_size'] + '_ad'
    elif(gd.variables['gas_flag'] == '2'):
        gd.mainDir =  'gadg_' + gd.variables['n'] + '_' + gd.variables['box_size'] + '_ch'
        
    if(not os.path.exists(gd.mainDir)):
        os.makedirs(gd.mainDir)
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################


if __name__=="__main__":
    go()


##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################   
