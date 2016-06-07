import os
import sys

json_path=os.path.abspath(sys.argv[1])
if not os.path.isfile(json_path):
    print "ERROR: Please provide path to a valid config.json file..."
    print sys.argv[1]
    exit(1)
    
import json
with open(json_path) as json_file:
    parameter = json.load(json_file)

    
if not parameter["BASE_DIR"][-1:] == '/':
    print "ERROR (BASE_DIR): Directories in config.json should end with '/'"
    print parameter["BASE_DIR"]
    exit(1)
if not parameter["OUTPUT_DIR"][-1:] == '/':
    print "ERROR (OUTPUT_DIR): Directories in config.json should end with '/'"
    print parameter["OUTPUT_DIR"]
    exit(1)
if not parameter["SAVE_DIR"][-1:] == '/':
    print "ERROR (SAVE_DIR): Directories in config.json should end with '/'"
    print parameter["SAVE_DIR"]
    exit(1)
if not parameter["SOURCE_DIR"][-1:] == '/':
    print "ERROR (SOURCE_DIR): Directories in config.json should end with '/'"
    print parameter["SOURCE_DIR"]
    exit(1)
if not parameter["kallisto"]["TCC_output"][-1:] == '/':
    print "ERROR (TCC_output): Directories in config.json should end with '/'"
    print parameter["kallisto"]["TCC_output"]
    exit(1)
    
    
JSON_ERR=0          
for i in range(len(parameter["read_filenames"])):
    if not os.path.isfile(str(parameter["BASE_DIR"])+str(parameter["read_filenames"][i])):
        JSON_ERR=1
        print "config.json error:"+str(parameter["BASE_DIR"])+str(parameter["read_filenames"][i])     
    

os.chdir(parameter["SOURCE_DIR"])
try:        
    if not os.path.isfile(str(parameter["SAVE_DIR"])+"barcodes.dat"): 
        print "ERROR:"+str(parameter["SAVE_DIR"])+"barcodes.dat not found"; exit(1)
    os.system('python error_correct_and_split.py '+json_path)
    if not os.path.isfile(str(parameter["OUTPUT_DIR"])+"umi_read_list.txt"): 
        print "ERROR:"+str(parameter["OUTPUT_DIR"])+"umi_read_list.txt not found"; exit(1)
    os.system('python compute_TCCs.py '+json_path)
    if not os.path.isfile(parameter["kallisto"]["TCC_output"]+"matrix.tsv"): 
        print "ERROR:"+parameter["kallisto"]["TCC_output"]+"matrix.tsv not found"; exit(1)    
    os.system('python prep_TCC_matrix.py '+json_path)
except:
    print "ERROR."