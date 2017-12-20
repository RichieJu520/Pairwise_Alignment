# Emboss_Needle_batch.py

@author: Feng Ju
@email: richieju520@gmail.com
The script was written and tested in python 2.7
"""

#print 'load module for use emboss and python on Euler'
import subprocess
#subprocess.call('module load gcc/4.8.2 gdc python/2.7.11',shell=True)
#subprocess.call('module load gcc/4.8.2 gdc emboss/6.5.7',shell=True)

from Bio.Emboss.Applications import NeedleCommandline
from Bio import SeqIO

def EmbossNeedle(seq_fname1,seq_fname2):
    needle_fname = "Emboss.Needle.txt"
    needle_cli = NeedleCommandline(asequence=seq_fname1, \
                                   bsequence=seq_fname2, \
                                   gapopen=10, \
                                   gapextend=0.5, \
                                   outfile=needle_fname)

    """This generates the needle file"""
    needle_cli() 
    """That parses the needle file, aln[0] and aln[1] contain the aligned
    first and second sequence in the usual format (e.g. - for a gap)"""
    #aln = AlignIO.read(needle_fname, "emboss")
    
    needle = open("Emboss.Needle.txt",'r')

    for line in needle:
        if line.find("Identity") > 0:
           num_i = line[line.find(":") + 1:line.find("(")].strip()
           identity = line[line.find("(") + 1:line.find("%")]
        if line.find("Similarity") > 0:
           num_s = line[line.find(":") + 1:line.find("(")].strip()
           similarity = line[line.find("(") + 1:line.find("%")]

    #print 'Identity:', identity, num_i
    #print 'Similarity:', similarity, num_s
    needle.close()
    return num_i, identity, num_s, similarity

Fasta1 = "ARGs-3389s-finalized.faa"
Fasta2 = "ARGs-combined-db.faa"
blast  = "map.blastp"

a, b = {}, {}
for rec in SeqIO.parse(Fasta1,'fasta'):
    a[str(rec.id)] = str(rec.seq)

for rec in SeqIO.parse(Fasta2,'fasta'):
    b[str(rec.id)] = str(rec.seq)
    
print len(a),' seqs in the 1st fasta'
print len(b),' seqs in the 1st fasta'

# Read map.blastp and make alignments
i=0
f=open(blast+'-Emboss-Needle.identity.csv','w')
f.write(','.join(['query','subject','identity', 'similarity', 'i-match','s-match','identity_l', 'similarity_l', 'i-match_l','s-match_l'])+'\n')

m, n = 0, 0
for line in open(blast,'r'):
    f1 = open('seq1.fa','w')
    f2 = open('seq2.fa','w')
    f3 = open('seq2l.fa','w')
    
    i+=1
    lis = line.rstrip().split('\t')


    query = lis[0]
    subject = lis[1]
    s_start = min(int(lis[-3]),int(lis[-4]))
    s_end =   max(int(lis[-3]),int(lis[-4]))

    try:
        f1.write('>'+lis[0]+'\n'+a[lis[0]].replace('*',''))
        f2.write('>'+lis[1]+'\n'+b[lis[1]].replace('*',''))
        f3.write('>'+lis[1]+'\n'+b[lis[1]].replace('*','')[(s_start-1):s_end])
        
        f1.close()
        f2.close()
        f3.close()
        m+=1
        
    except KeyError:
        n+=1
        continue
    
    num_i, identity, num_s, similarity = EmbossNeedle('seq1.fa','seq2.fa')
    num_i_l, identity_l, num_s_l, similarity_l = EmbossNeedle('seq1.fa','seq2l.fa')
    
    f.write(','.join([query, subject, identity, similarity, num_i, num_s,identity_l, similarity_l, num_i_l, num_s_l])+'\n')

print m, 'pair of sequences calculated!'
print n, 'hits ignored!' 
print 'DONE!'
    
    
   
