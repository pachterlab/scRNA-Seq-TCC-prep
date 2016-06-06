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

os.chdir(parameter["SOURCE_DIR"])
try:        
    os.system('python error_correct_and_split.py '+json_path)
    os.system('python compute_TCCs.py '+json_path)
    os.system('python prep_TCC_matrix.py '+json_path)
except:
    print "ERROR."