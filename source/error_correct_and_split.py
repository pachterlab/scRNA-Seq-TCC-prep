
# coding: utf-8

import os
import sys
json_path=os.path.abspath(sys.argv[1])
if not os.path.isfile(json_path):
    print "ERROR: Please provide path to a valid config.json file..."
    print sys.argv[1]
    exit(1)

    
workpath, jsonfile = os.path.split(json_path)
os.chdir(workpath)

# In[1]:
import numpy as np
from itertools import islice
from collections import Counter
import matplotlib.pyplot as plt
import gzip, time, random
from multiprocessing import Pool
import pickle


# ### Read config.json file

# In[2]:

import json
with open("config.json") as json_file:
    parameter = json.load(json_file)


print "READ FILES:\n"
read_dirs=[]
for i in range(len(parameter["read_filenames"])):
    read_dirs+=[str(parameter["BASE_DIR"])+str(parameter["read_filenames"][i])]
    print read_dirs[i]
    
    
    
random.seed()
BARCODE_LENGTH=parameter['BARCODE_LENGTH']
output_dir = parameter['OUTPUT_DIR']
NUM_THREADS = parameter['NUM_THREADS']
#temporary file to extract all reads
all_reads_file = str(parameter["BASE_DIR"])+'all_reads.fastq'


# ### Load barcodes

# In[3]:

def encoding_map(ch):
    if ch=='A':return 0
    if ch=='G':return 1
    if ch=='C':return 2
    if ch=='T':return 3
    if ch=='N':return random.randint(0,3)

decoding_lst = ['A', 'G', 'C', 'T']

def encode(k):
    code = 0
    for ch in k:
        code *= 4
        code += encoding_map(ch)
    return code

def decode(code):
    ret = ''
    for _ in range(14):
        index = code & 3
        code >>= 2
        ret = decoding_lst[index] + ret
    return ret


#LOAD barcodes
save_dir=str(parameter["SAVE_DIR"])

print "Loading Barcodes..."
t0 = time.time()
with open(save_dir+"barcodes.dat", 'rb') as f:
    barcodes=pickle.load(f)
with open(save_dir+"codewords.dat", 'rb') as f:
    codewords=pickle.load(f)
with open(save_dir+"brc_idx_to_correct.dat", 'rb') as f:
    brc_idx_to_correct= pickle.load(f)
t1 = time.time()
print t1-t0, "sec"


# ### Error-correct barcodes

# In[4]:

def merge_barcodes(id):
    if id in brc_idx_to_correct:
        s=set(mutations(decode(codewords[id]),1))
        pos=[]
        for idx, barcode in enumerate(barcodes):
            if barcode in s:
                pos+=[idx]
        return pos
    else:
        s=codewords[id]
        return [idx for idx, barcode in enumerate(barcodes) if barcode == s]
    
    
import itertools
def mutations(word, hamming_distance, charset='ATCG'):
    for indices in itertools.combinations(range(len(word)), hamming_distance):
        for replacements in itertools.product(charset, repeat=hamming_distance):
            mutation = list(word)
            for index, replacement in zip(indices, replacements):
                mutation[index] = replacement
            yield encode("".join(mutation))



print "Merging barcodes... NUM_THREADS =",NUM_THREADS 
p = Pool(NUM_THREADS)
t0 = time.time()
ret_vec=p.map(merge_barcodes, range(len(codewords)))
t1 = time.time()
print t1-t0, "sec"
p.close(); p.join()


reads_per_barcode=[]
for i in range(len(codewords)):
    reads_per_barcode+=[len(ret_vec[i])]
print "Reads in Barcodes:",sum(reads_per_barcode)


# ### Output single-cell files

# In[5]:

#create output directory 
import os
if not os.path.isdir(output_dir):
    try:
        os.mkdir(output_dir)
    except OSError as e:
        print "OSError({0}): {1}".format(e.errno, e.strerror)
    


# In[6]:

#concatenate all .gz read files

command = "cat "
for files in [read_dirs[0],read_dirs[1],read_dirs[2],read_dirs[3],read_dirs[4],read_dirs[5],read_dirs[6],read_dirs[7]]:  
    command+=files+' '
command+="> "+all_reads_file+".gz"
print "cat..."
os.system(command)


# In[7]:

# temporarilly unzip all reads

t0=time.time()
print "gunzip..."

os.system("gunzip "+all_reads_file+".gz")

t1=time.time()
print t1-t0, "sec"


# In[8]:

# create line_offset list

f = open(all_reads_file)

t0=time.time()
print "line_offset..."
line_offset = []
offset = 0
for line in f:
    line_offset.append(offset)
    offset += len(line)
    
f.close()
t1=time.time()
print t1-t0, "sec"  


NUM_OF_LINES=len(line_offset)
print "number of reads in dataset =",NUM_OF_LINES/8


# In[9]:

# Split single-cell files and umis

f = open(all_reads_file)
t0=time.time()
for cell in range(len(codewords)):
    filename = "cell_"+str(cell).zfill(4)+'_'+decode(codewords[cell])
    print "writing " + filename +"..."
    output_umis=""
    output_fastq=""
    for i in ret_vec[cell]:
        f.seek(line_offset[i*8])
        output_fastq+=f.readline()
        output_fastq+=f.readline()
        output_fastq+=f.readline()
        output_fastq+=f.readline()

        f.seek(line_offset[5+i*8])
        output_umis+=f.readline()
    
    with open(output_dir+filename+".umi", 'wb') as umi:
        umi.write(output_umis)
    with open(output_dir+filename+".fastq", 'wb') as reads:
        reads.write(output_fastq)

f.close()
t1=time.time()
print t1-t0, "sec" 


# In[10]:

#remove temp all_reads file 
os.system("rm "+all_reads_file)


#compress output files 
from os import listdir
from os.path import isfile, join
fastqfiles = [output_dir+f for f in listdir(output_dir) if isfile(join(output_dir, f)) and f[-6:]==".fastq"]

def gzip_fastqs(filepath):
    if filepath[-6:]==".fastq":
        os.system("gzip "+ filepath)

print "gzip..."

p=Pool(8)
t0 = time.time()
p.map(gzip_fastqs, fastqfiles)
t1 = time.time()
print t1-t0, "sec"
p.close(); p.join()


# In[11]:

#create batch file: singlecell_umi_read_list.txt
fastqfiles = [output_dir+f for f in sorted(listdir(output_dir)) if isfile(join(output_dir, f)) and f[-9:]==".fastq.gz"]
umifiles = [output_dir+f for f in sorted(listdir(output_dir)) if isfile(join(output_dir, f)) and f[-4:]==".umi"]
cell_ids = [f[:24] for f in sorted(listdir(output_dir)) if isfile(join(output_dir, f)) and f[-4:]==".umi"]

out_data=''
for i in range(len(cell_ids)):
    out_data+=cell_ids[i]+'\t'+umifiles[i]+'\t'+fastqfiles[i]+'\n'

with open(str(parameter["BASE_DIR"])+"umi_read_list.txt", 'wb') as f:
                   f.write(out_data)
    
print "DONE"    


# In[ ]:




# In[ ]:



