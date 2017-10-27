###### v0.2 ######
# coding: utf-8

import os
import sys
json_path=os.path.abspath(sys.argv[1])
if not os.path.isfile(json_path):
    print("ERROR: Please provide path to a valid config.json file...")
    print(sys.argv[1])
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

from os import listdir
from os.path import isfile, join
barcode_filenames = [f for f in sorted(listdir(str(parameter["BASE_DIR"]))) if isfile(join(str(parameter["BASE_DIR"]), f)) and f[:7]=="read-I1" and f[11:19] in parameter["sample_idx"]] 
read_filenames = ['read-RA'+f[7:] for f in barcode_filenames]

print("READ FILES:\n")
read_dirs=[]
for i in range(len(read_filenames)):
    read_dirs+=[str(parameter["BASE_DIR"])+read_filenames[i]]
    print(read_dirs[i])
    
    
    
random.seed()
BARCODE_LENGTH=parameter['BARCODE_LENGTH']
output_dir = parameter['OUTPUT_DIR']
NUM_THREADS = parameter['NUM_THREADS']
#temporary file to extract all reads
all_reads_file = str(parameter["SAVE_DIR"])+'all_reads.fastq'


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
    for _ in range(BARCODE_LENGTH):
        index = code & 3
        code >>= 2
        ret = decoding_lst[index] + ret
    return ret


#LOAD barcodes
save_dir=str(parameter["SAVE_DIR"])

print("Loading Barcodes...")
t0 = time.time()
with open(save_dir+"barcodes.dat", 'rb') as f:
    barcodes=pickle.load(f)
with open(save_dir+"codewords.dat", 'rb') as f:
    codewords=pickle.load(f)
with open(save_dir+"brc_idx_to_correct.dat", 'rb') as f:
    brc_idx_to_correct= pickle.load(f)
t1 = time.time()
print(t1-t0, "sec")


# ### Error-correct barcodes

# In[4]:
chunksize=1+int(len(barcodes)/NUM_THREADS)

cw={}
for id in range(len(codewords)):
    cw[codewords[id]] = id

barcode_split=[]
for i in range(0, len(barcodes), chunksize):        
    barcode_split+=[[i,barcodes[i:i+chunksize]]]

from itertools import chain, combinations, product
def hamming_circle(s, n, alphabet='ATCG'):
    """Generate strings over alphabet whose Hamming distance from s is
    exactly n.
    """
    for positions in combinations(list(range(len(s))), n):
        for replacements in product(list(range(len(alphabet) - 1)), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield encode(''.join(cousin))
            
def merge_barcodes(barcs):
    offset=barcs[0]
    barcs=barcs[1]
    retvec=[]
    for id in range(len(codewords)):
        retvec+=[[]]
    for idx, barcode in enumerate(barcs):
        if barcode in codeword_set: retvec[cw[barcode]] +=[idx+offset]
        else:
            neighbors = hamming_circle(decode(barcode),1)
            for neighbor in neighbors:
                if neighbor in brc_to_correct: retvec[cw[neighbor]] +=[idx+offset]; break;
    return retvec
            
         
print("Merging barcodes...") 
codeword_set = set(codewords)
codeword_list = list(codewords)
brc_to_correct=set(codewords[brc_idx_to_correct])

t0 = time.time()
p = Pool(NUM_THREADS)
t0 = time.time()
ret_threads=p.map(merge_barcodes, barcode_split)
p.close(); p.join()


ret_vec=[]
for id in range(len(codewords)):
    idx_list=[]
    for t in range(len(ret_threads)):
        idx_list+=ret_threads[t][id]        
    ret_vec+=[idx_list]

#ret_vec=[]
#for id in range(len(codewords)):
#    idx_list=[item for sublist in ret_threads[:][id] for item in sublist]       
#    ret_vec+=[idx_list]

    
t1 = time.time()
print(t1-t0, "sec")

reads_per_barcode=[]
for i in range(len(codewords)):
    reads_per_barcode+=[len(ret_vec[i])]
NUM_OF_READS_in_CELL_BARCODES = sum(reads_per_barcode)
print("NUM_OF_READS_in_CELL_BARCODES (after error-correct):",NUM_OF_READS_in_CELL_BARCODES) 

# ### Output single-cell files

# In[5]:

#create output directory 
import os
if not os.path.isdir(output_dir):
    try:
        os.mkdir(output_dir)
    except OSError as e:
        print("OSError({0}): {1}".format(e.errno, e.strerror))
    


# In[6]:

#concatenate all .gz read files

command = "cat "
for files in read_dirs:  
    command+=files+' '
command+="> "+all_reads_file+".gz"
print("cat...")
os.system(command)


# In[7]:

# temporarilly unzip all reads

t0=time.time()
print("gunzip...")

os.system("gunzip -f "+all_reads_file+".gz")

t1=time.time()
print(t1-t0, "sec")


# In[8]:

# create line_offset list

f = open(all_reads_file)

t0=time.time()
print("line_offset...")
line_offset = []
offset = 0
for line in f:
    line_offset.append(offset)
    offset += len(line)
    
f.close()
t1=time.time()
print(t1-t0, "sec")  


NUM_OF_LINES=len(line_offset)
print("number of reads in dataset =",NUM_OF_LINES/8)


# In[9]:

# Split single-cell files and umis

f = open(all_reads_file)
t0=time.time()
for cell in range(len(codewords)):
    filename = "cell_"+str(cell).zfill(4)+'_'+decode(codewords[cell])
    print("writing " + filename +"...")
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
print(t1-t0, "sec") 


# In[10]:

#remove temp all_reads file 
os.system("rm "+all_reads_file)


#compress output files 
from os import listdir
from os.path import isfile, join
fastqfiles = [output_dir+f for f in listdir(output_dir) if isfile(join(output_dir, f)) and f[-6:]==".fastq"]

def gzip_fastqs(filepath):
    if filepath[-6:]==".fastq":
        os.system("gzip -f "+ filepath)

print("gzip...")

p=Pool(8)
t0 = time.time()
p.map(gzip_fastqs, fastqfiles)
t1 = time.time()
print(t1-t0, "sec")
p.close(); p.join()


# In[11]:

#create batch file: singlecell_umi_read_list.txt
fastqfiles = [output_dir+f for f in sorted(listdir(output_dir)) if isfile(join(output_dir, f)) and f[-9:]==".fastq.gz"]
umifiles = [output_dir+f for f in sorted(listdir(output_dir)) if isfile(join(output_dir, f)) and f[-4:]==".umi"]
cell_ids = [f[:24] for f in sorted(listdir(output_dir)) if isfile(join(output_dir, f)) and f[-4:]==".umi"]

out_data=''
for i in range(len(cell_ids)):
    out_data+=cell_ids[i]+'\t'+umifiles[i]+'\t'+fastqfiles[i]+'\n'

with open(str(parameter["OUTPUT_DIR"])+"umi_read_list.txt", 'wb') as f:
                   f.write(out_data)

printer=""
printer+="NUM_OF_READS_in_CELL_BARCODES (after error-correct): %s\n" % NUM_OF_READS_in_CELL_BARCODES
printer+="NUMBER_OF_LINES in 'all_reads.fastq': %s\n" % NUM_OF_LINES  
with open(save_dir+"run.info", 'a') as f:
    f.write(printer)
print('\n')
print(printer)        
        
print("DONE")    

