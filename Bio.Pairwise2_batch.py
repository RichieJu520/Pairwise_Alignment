# -*- coding: utf-8 -*-
"""
@author: Feng Ju
@email: richieju520@gmail.com
The script was written and tested in python 2.7 and biopython 1.58
"""

from Bio import pairwise2
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo as matlist

def pairwise(asequence, bsequence):
    
    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -0.5
    
    alns = pairwise2.align.globalds(asequence, bsequence, matrix, gap_open, gap_extend)

    top_aln = alns[0]
    aln_a, aln_b, score, start, end = top_aln
    alignment = aln_a+'\n'+aln_b

    # calculate identity
    A=list(aln_a)
    B=list(aln_b)
    count, gaps = 0, 0

    len_A = len(aln_a.replace('-',''))
    len_B = len(aln_b.replace('-',''))

    for n in range(0, len(A)):
        if A[n]==B[n]:
            if B[n]!="-":
                count=count+1
            else:
                gaps=gaps+1
        else:
            continue
    length = end - start
    Identity = 100*count/float(len(A)-gaps)
    return Identity, score, length

# Read two fasta files
Fasta1 = "ARGs-3389s-finalized.faa"
Fasta2 = "ARGs-combined-db.faa"
blast  = "map.blastp"

a, b = {}, {}
for rec in SeqIO.parse(Fasta1,'fasta'):
    a[str(rec.id)] = str(rec.seq)

for rec in SeqIO.parse(Fasta2,'fasta'):
    b[str(rec.id)] = str(rec.seq)
    
print len(a),' seqs in the 1st fasta'
print len(b),' seqs in the 2nd fasta'

# Read map.blastp and make alignments
i=0
f=open(blast+'-Pairwise2-Needle-identity.csv','w')
f.write(','.join(['query','subject','identity','score','length','identity_l','score_l','length_l'])+'\n')

m, n = 0, 0
for line in open(blast,'r'):
    i+=1
    if i%1000==0:
        print i, 'pairs of seqs aligned!'
    lis = line.rstrip().split('\t')
    try:
        query = a[lis[0]].replace('*','')
        subject = b[lis[1]].replace('*','')
        s_start = min(int(lis[-3]),int(lis[-4]))
        s_end =   max(int(lis[-3]),int(lis[-4]))
        subject_l = b[lis[1]].replace('*','')[(s_start-1):s_end]
        Identity, score, length = pairwise(query, subject)
        Identity_l, score_l, length_l = pairwise(query, subject_l)
        m+=1
    except KeyError:
        n+=1
        continue
    f.write(','.join([lis[0],lis[1],str(Identity), str(score),str(length),str(Identity_l), str(score_l),str(length_l)])+'\n')
    
print m, 'pair of sequences calculated!'
print n, 'hits ignored!' 
print 'DONE'
    
                              

                    
