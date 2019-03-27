import os
import numpy
import scipy.sparse as sps
import scipy.io as scio
import sys
import time

start_time = time.time()

dir = next(os.walk('.'))[1]

for curr in dir:
    subdir = next(os.walk(curr))[1]

    for sub_curr in subdir:
        start_time_this = time.time()
        with open(curr + "/" + sub_curr + "/" + sub_curr.split("_STAR")[0] + ".cnts", 'w') as outfile, open("candidates_non_overlapping.txt", 'r') as RNA_list:
            
            curr_chr = RNA_list.readline().rstrip()
            while True:
                if curr_chr =='':
                    break
            
                curr_algn = scio.mmread(curr + "/" + sub_curr + "/" + curr_chr + ".mtx").tocsr()
            
                RNA_list.readline()
                RNA_list.readline()
                
                next_gene = RNA_list.readline()
                while True:
                    if next_gene == '\n':
                        break
                    gene_id, start, end = next_gene.rstrip().split('\t')
                    start = int(start)
                    end = int(end)
                    length = end-start +1
                    curr_gene= numpy.concatenate([curr_algn[nuc, (start-1):(end)].todense() for nuc in range(start-1, end)])
                    """if ('anti' in curr_chr) :
                        curr_gene = numpy.transpose(curr_gene)"""
                
                    if ('anti' not in curr_chr):
                        curr_counts = [numpy.sum(curr_gene[nuc+1, :]) for nuc in range(length-1)]
                        curr_cov = [numpy.sum(curr_gene[0:(nuc+2), (nuc+1):]) for nuc in range(length-1)]
                        curr_counts.append(0)
                        curr_cov.append(0)
                    else:
                        curr_counts = [numpy.sum(curr_gene[:, nuc-1]) for nuc in range(1, length)]
                        curr_cov = [numpy.sum(curr_gene[0:nuc, (nuc-1):length]) for nuc in range(1, length)]
                        curr_counts = numpy.concatenate(([0], curr_counts))
                        curr_cov = numpy.concatenate(([0], curr_cov))
                    
                    outfile.write(gene_id +"\n" + '\t'.join(["%d" % x for x in curr_counts]) + '\n' + '\t'.join(["%d" % x for x in curr_cov]) + '\n\n')
                    next_gene = RNA_list.readline()
                RNA_list.readline()
                curr_chr = RNA_list.readline().rstrip()
        print(curr + "/" + sub_curr + " processed in --- %s seconds ---" % (time.time() - start_time_this))

print("All processed in --- %s seconds ---" % (time.time() - start_time))