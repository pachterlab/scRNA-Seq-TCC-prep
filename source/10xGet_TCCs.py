import os
import sys

json_path=os.path.abspath(sys.argv[1])
if not os.path.isfile(json_path):
    print("ERROR: Please provide path to a valid config.json file...")
    print(sys.argv[1])
    exit(1)
    
import json
with open(json_path) as json_file:
    parameter = json.load(json_file)
    

if not parameter["OUTPUT_DIR"][-1:] == '/':
    print("ERROR (OUTPUT_DIR): Directories in config.json should end with '/'")
    print(parameter["OUTPUT_DIR"])
    exit(1)
if not parameter["SAVE_DIR"][-1:] == '/':
    print("ERROR (OUTPUT_DIR): Directories in config.json should end with '/'")
    print(parameter["SAVE_DIR"])
    exit(1)
if not parameter["SOURCE_DIR"][-1:] == '/':
    print("ERROR (OUTPUT_DIR): Directories in config.json should end with '/'")
    print(parameter["SOURCE_DIR"])
    exit(1)
if not parameter["kallisto"]["TCC_output"][-1:] == '/':
    print("ERROR (OUTPUT_DIR): Directories in config.json should end with '/'")
    print(parameter["kallisto"]["TCC_output"])
    exit(1)
for path in parameter["FASTQ_DIRS"]:
    if not path[-1:] == '/':
        print("ERROR (OUTPUT_DIR): Directories in config.json should end with '/'")
        print(path)
        exit(1)
    if not os.path.isdir(path):
        print("config.json error: FASTQ_DIRS contain invalid paths")
        exit(1)
        
if not os.path.isdir(parameter["SOURCE_DIR"]):
    print("config.json error: SOURCE_DIR is not a valid path")
    exit(1)    
    
try:
    os.system(parameter["kallisto"]["binary"]+" version")
except:
    print("ERROR with kallisto binary "+parameter["kallisto"]["binary"])
    exit(1)

    
########################################## check json file parameters
from os import listdir
from os.path import isfile, join
import re


FASTQ_DIRS=parameter["FASTQ_DIRS"]
SAMPLE_NAMES=parameter['sample_names']
EXP_CELLS=parameter['EXP_CELLS']

if len(EXP_CELLS)!=len(SAMPLE_NAMES):
    print("ERROR: Expected cells must be provided for every sample")
    
try:
    tmp=SAMPLE_NAMES[0]
except TypeError:
    print("Sample names must be provided as a list (even when there is only one sample, e.g., 'sample1' ---> ['sample1'] ). ")

try:
    tmp=EXP_CELLS[0]
except TypeError:
    print("EXP_CELLS must be provided as a list (even when there is only one sample, e.g., 3000 ---> [3000] ). ")

print("\n\nWill process samples:")
for i in range(len(SAMPLE_NAMES)):
    print(' '+SAMPLE_NAMES[i]+" with ~%d cells"% EXP_CELLS[i])
    
print("\nWill look for the corresponding fastqs in the directories:")    
for d in FASTQ_DIRS:
    print(' '+d)
    
##################################
try:    
    barcode_filenames_per_sample = []
    barcode_dirs_per_sample = []
    read_filenames_per_sample = []
    read_dirs_per_sample = []
    for i_S in range(len(SAMPLE_NAMES)):
        barcode_filenames = []
        barcode_dirs = []
        read_filenames = []
        read_dirs = []
        for i_F in range(len(FASTQ_DIRS)):
            for f in sorted(listdir(str(FASTQ_DIRS[i_F]))):
                
                if (isfile(join(str(FASTQ_DIRS[i_F]), f)) and f.startswith(str(SAMPLE_NAMES[i_S])) and 
                    f[-8:]=='fastq.gz' and f.find('R1')>0):

                    barcode_filenames.append(f)
                    barcode_dirs.append(str(FASTQ_DIRS[i_F])+f)

                if (isfile(join(str(FASTQ_DIRS[i_F]), f)) and f.startswith(str(SAMPLE_NAMES[i_S])) and 
                    f[-8:]=='fastq.gz' and f.find('R2')>0):

                    read_filenames.append(f)
                    read_dirs.append(str(FASTQ_DIRS[i_F])+f)
                    
#                pattern = r"^"+str(SAMPLE_NAMES[i_S])+"_S\d+_L00\d_([IR]\d)_001.fastq.gz"
#                match = re.match(pattern, f)
#                if match:
#                    barcode_name = match.group(1)
#                if isfile(join(str(FASTQ_DIRS[i_F]), f)) and f.startswith(str(SAMPLE_NAMES[i_S])) and barcode_name == 'R1':
#                    barcode_filenames.append(f)
#                    barcode_dirs.append(str(FASTQ_DIRS[i_F])+f)
#                if isfile(join(str(FASTQ_DIRS[i_F]), f)) and f.startswith(str(SAMPLE_NAMES[i_S])) and barcode_name == 'R2':
#                    read_filenames.append(f)
#                    read_dirs.append(str(FASTQ_DIRS[i_F])+f)
        barcode_filenames_per_sample += [barcode_filenames]
        barcode_dirs_per_sample += [barcode_dirs]
        read_filenames_per_sample += [read_filenames]
        read_dirs_per_sample +=[read_dirs]
    print("\n----------------------------------------------------------------------------------------")               
    print("BARCODE FILES(R1):")
    for i in range(len(barcode_dirs_per_sample)):
        print("\nsample: "+SAMPLE_NAMES[i])
        for j in range(len(barcode_dirs_per_sample[i])):
            print(" "+barcode_dirs_per_sample[i][j])    
    print("\n----------------------------------------------------------------------------------------")
    print("READ FILES(R2):")
    for i in range(len(read_dirs_per_sample)):
        print("\nsample: "+SAMPLE_NAMES[i])
        for j in range(len(read_dirs_per_sample[i])):
            print(" "+read_dirs_per_sample[i][j])    

    print("\n----------------------------------------------------------------------------------------")   
except:
    print("ERROR loading R1 and/or R2 files for samples: "+str(SAMPLE_NAMES))
    
os.chdir(parameter["SOURCE_DIR"])
skip=''
try:
    if os.path.isdir(str(parameter["SAVE_DIR"])):
        if os.path.isfile(str(parameter["SAVE_DIR"])+"notebook.cookie"):
            f= open(str(parameter["SAVE_DIR"])+"notebook.cookie", 'r')
            skip=f.readline()
            f.close()
    if skip==str(SAMPLE_NAMES):
        print("It seems that the 10xDetect_cell_barcodes notebook has already been executed for samples "+skip)
        print("Will SKIP get_cell_barcodes script (you can remove the file '"+str(parameter["SAVE_DIR"])+"notebook.cookie' to re-run this step...)")
    else:
        os.system('python get_cell_barcodes_v2.py '+json_path)
    if not os.path.isfile(str(parameter["SAVE_DIR"])+SAMPLE_NAMES[0]+"_barcodes.dat"): 
        print("ERROR:"+str(parameter["SAVE_DIR"])+SAMPLE_NAMES[0]+"_barcodes.dat not found"); exit(1)
    os.system('python error_correct_and_split_v2.py '+json_path)
    if not os.path.isfile(str(parameter["OUTPUT_DIR"])+"umi_read_list.txt"): 
        print("ERROR:"+str(parameter["OUTPUT_DIR"])+"umi_read_list.txt not found"); exit(1)
    os.system('python compute_TCCs.py '+json_path)
except:
    print("ERROR.") 

    
    
    
    
    