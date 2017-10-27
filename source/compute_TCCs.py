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

kallisto_cmd= parameter["kallisto"]["binary"]+" pseudo -i "+parameter["kallisto"]["index"]+" -o "+parameter["kallisto"]["TCC_output"]+" --umi -b "+parameter["OUTPUT_DIR"]+"umi_read_list.txt"+" -t "+str(parameter["NUM_THREADS"])

print("Running kallisto pseudo:")
print(kallisto_cmd)

printer="KALLISTO CMD: "+kallisto_cmd+'\n'

save_dir=str(parameter["SAVE_DIR"])
with open(save_dir+"run.info", 'a') as f:
    f.write(printer)

os.system(kallisto_cmd)
print("DONE.")