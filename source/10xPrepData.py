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

    
if not parameter["BASE_DIR"][-1:] == '/':
    print("ERROR (BASE_DIR): Directories in config.json should end with '/'")
    print(parameter["BASE_DIR"])
    exit(1)
if not parameter["OUTPUT_DIR"][-1:] == '/':
    print("ERROR (OUTPUT_DIR): Directories in config.json should end with '/'")
    print(parameter["OUTPUT_DIR"])
    exit(1)
if not parameter["SAVE_DIR"][-1:] == '/':
    print("ERROR (SAVE_DIR): Directories in config.json should end with '/'")
    print(parameter["SAVE_DIR"])
    exit(1)
if not parameter["SOURCE_DIR"][-1:] == '/':
    print("ERROR (SOURCE_DIR): Directories in config.json should end with '/'")
    print(parameter["SOURCE_DIR"])
    exit(1)
if not parameter["kallisto"]["TCC_output"][-1:] == '/':
    print("ERROR (TCC_output): Directories in config.json should end with '/'")
    print(parameter["kallisto"]["TCC_output"])
    exit(1)
    
from os import listdir
from os.path import isfile, join
barcode_filenames = [f for f in sorted(listdir(str(parameter["BASE_DIR"]))) if isfile(join(str(parameter["BASE_DIR"]), f)) and f[:7]=="read-I1" and f[11:19] in parameter["sample_idx"]]
read_filenames = ['read-RA'+f[7:] for f in barcode_filenames]

if len(barcode_filenames)==0:
    JSON_ERR=1
    print("ERROR: no barcode files (read-I1) were found in" + str(parameter["BASE_DIR"]))

JSON_ERR=0          
for i in range(len(read_filenames)):
    if not os.path.isfile(str(parameter["BASE_DIR"])+read_filenames[i]):
        JSON_ERR=1
        print("config.json error:"+str(parameter["BASE_DIR"])+read_filenames[i])       

os.chdir(parameter["SOURCE_DIR"])
try:        
    if not os.path.isfile(str(parameter["SAVE_DIR"])+"barcodes.dat"): 
        print("ERROR:"+str(parameter["SAVE_DIR"])+"barcodes.dat not found"); exit(1)
    os.system('python error_correct_and_split.py '+json_path)
    if not os.path.isfile(str(parameter["OUTPUT_DIR"])+"umi_read_list.txt"): 
        print("ERROR:"+str(parameter["OUTPUT_DIR"])+"umi_read_list.txt not found"); exit(1)
    os.system('python compute_TCCs.py '+json_path)
    if not os.path.isfile(parameter["kallisto"]["TCC_output"]+"matrix.tsv"): 
        print("ERROR:"+parameter["kallisto"]["TCC_output"]+"matrix.tsv not found"); exit(1)    
    os.system('python prep_TCC_matrix.py '+json_path)
except:
    print("ERROR.")