# Emboss_Stretcher_batch.py

#print 'load module for use emboss and python on Euler'
import subprocess
#subprocess.call('module load gcc/4.8.2 gdc python/2.7.11',shell=True)
#subprocess.call('module load gcc/4.8.2 gdc emboss/6.5.7',shell=True)

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
f=open(blast+'-Emboss-Stretcher.identity.csv','w')
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

    subprocess.call('needle seq1.fa seq2.fa Emboss1.Needle.txt',shell=True)  ### use the Emboss's Stretcher installed in the linux system
    subprocess.call('needle seq1.fa seq2l.fa Emboss2.Needle.txt',shell=True) ### use the Emboss's Stretcher installed in the linux system

    needle1 = open("Emboss.Needle1.txt",'r')
    for line in needle1:
        if line.find("Identity") > 0:
           num_i = line[line.find(":") + 1:line.find("(")].strip()
           identity = line[line.find("(") + 1:line.find("%")]
        if line.find("Similarity") > 0:
           num_s = line[line.find(":") + 1:line.find("(")].strip()
           similarity = line[line.find("(") + 1:line.find("%")]

    needle2 = open("Emboss.Needle2.txt",'r')
    for line in needle2:
        if line.find("Identity") > 0:
           num_i_l = line[line.find(":") + 1:line.find("(")].strip()
           identity_l = line[line.find("(") + 1:line.find("%")]
        if line.find("Similarity") > 0:
           num_s_l = line[line.find(":") + 1:line.find("(")].strip()
           similarity_l = line[line.find("(") + 1:line.find("%")]
           
    needle1.close()
    needle2.close()
    f.write(','.join([query, subject, identity, similarity, num_i, num_s,identity_l, similarity_l, num_i_l, num_s_l])+'\n')

print 'DONE'
