#!/usr/bin/env python3
# _*_ coding: utf-8 _*_
 
# @File    :   homolog.py
# @Version :   3.0 portable
# @Author  :   NieYuqi
# @Email   :   nieyuqi@163.com
# @Time    :   2021/05/31 16:05:19

#Description:
    # This script use the hmmscan to find the similar gene to the query gene, and use the blastp to polish them.
#Update:
    # 2021/06/01 15:36:00
        # v2.0 add a mechanism to test the query gene in the inital input file if similar.
    # 2021/06/01 16:10:00
        # v3.0 add a mechanism to dislodge the sequences whose ranking is non-top in Hmmscan tblout.

import datetime
start_time = datetime.datetime.now()
print("{0:=^40}".format(' Start '))

import os  #os.system()
import sys #sys.argv[1]

# 参数设置
threadsnum = str(4) #blastp hmmscan mafft cpu核心数 
sourcepath = './source/'
referencedb = sourcepath+'Arabidopsis_thaliana.fa'

# create dir
temppath = './temp/'
outpath = './out/'
if os.path.exists("temp"):
    os.rename("temp","temp_old_"+str(start_time))
    print(" Warning ".center(40,"-"))
    print("[Warning]:The file name has been changed: 'temp'"+" to "+"'temp_old_"+str(start_time).replace(" ","_")+"'.")
os.system('mkdir temp')
if os.path.exists("out"):
    os.rename("out","out_old_"+str(start_time))
    print(" Warning ".center(40,"-"))
    print("[Warning]:The file name has been changed: 'out'"+" to "+"'out_old_"+str(start_time).replace(" ","_")+"'.")
os.system('mkdir out')

# # extract sequence from fasta
# def extract_fa(id,fa,out):
# 	dic = []
# 	for line in fa:
# 		line = line.strip()
# 		if line.startswith(">"):
# 			seqid = line.split()[0][1:]
# 			dic[seqid] = ''
# 		else:
# 			dic[seqid] = dic.get(seqid,'') + line
# 	for id,seq in dic.items():

# extract sequence from fasta flie
def extract_seq(id_list,fasta,outfile): #三个参数分别为：序列号(列表)，输入fasta，输出fasta
    if id_list != None:
        with open(fasta,'r') as fa:
            dic = {}
            for line in fa:
                line = line.strip()
                if line.startswith('>'):
                    # seq_id = line.strip().split("\t")[0]
                    seq_id = line.split()[0]
                    seq_id = seq_id[1:]
                    dic[seq_id] = ''
                else:
                    dic[seq_id] += line
                #    dic[seq_id] += line
        with open(outfile,'w') as out:
            for i in id_list:
                seq=dic[i]
                out.write('>'+i+'\n'+seq+'\n')


# align
def align(alignin,alignout):
    global threadsnum
    os.system('./source/mafft --auto --thread '+threadsnum+' '+alignin+'>'+alignout)

# hmm prepare
def hmmpre(fasta,hmmfile,hmmname):
    os.system('./source/hmmbuild  -n '+hmmname+" "+hmmfile+" "+fasta)
    os.system('./source/hmmpress '+hmmfile)

# hmmscan
def hmm(hmmfile,hmmout,hmmtblout,fasta):
    global threadsnum
    os.system('./source/hmmscan --cpu '+threadsnum+' -o '+hmmout+' --tblout '+hmmtblout+' '+hmmfile+' '+fasta )

# analyse hmmscan
def analyse(hmmtblout): # return a lists ranked by the best 1 domin score(candidate query)
    candidate = []
    with open(hmmtblout,'r') as tblout:
        for line in tblout:
            line = line.strip()
            if line.startswith('#'):
                continue
            else:
                line=line.split()
                ls=[]
                id = line[2]
                score = line[8] # best 1 domin score
                ls.append(id)
                ls.append(score)
                candidate.append(ls)
    candidate.sort(key=lambda x:eval(x[1]),reverse = True)
    candi = []
    for i in candidate:
        candi.append(i[0])
    return candi,candidate

# blastp
def blast(query,db,blastout,targetnum):
    global threadsnum
    targetnum = str(targetnum)
    os.system('./source/blastp -query '+query+' -db '+db+' -out '+blastout+' -outfmt 6 -max_target_seqs '+targetnum+'  -num_threads ' + threadsnum)

# trace blastp out
def trace(blastout,ids,traceout):
    with open(blastout,'r') as txt:
        dic = {}
        for line in txt:
            line = line.strip()
            line = line.split()
            query = line[0]
            target = line[1]
            dic[query] = dic.get(query,[])
            dic[query].append(target)
    with open(traceout,'w') as out:
        for query,target in dic.items():
            iftarget = True
            for id in target:
                if iftarget == False:
                    break
                else:
                    if id in ids:
                        continue
                    elif id not in ids:
                        iftarget = False
            if iftarget == True:
                out.write(query+'\n')

# test the if the input locus are homology
def testhomo(idlist,genename):
    query = temppath+genename+".fa"
    blastout = temppath+genename+'.test.blast'
    targetnum = len(idlist)
    blast(query,referencedb,blastout,targetnum)
    traceout = temppath+genename+".trace.txt"
    trace(blastout,idlist,traceout)
    with open(traceout,'r') as txt:
        ls = []
        for line in txt:
            line = line.strip()
            if line!='':
                ls.append(line)
        if len(ls) == len(idlist):
            testhomopass = True
        else:
            testhomopass = False
        for i in ls:
            if i not in idlist:
                testhomopass = False
                break
    return testhomopass

# dislodge the sequence
def dislodge(traceout,candidatels2,dislodgeout):
    with open(traceout,'r') as read:
        ls = []
        for line in read:
            line = line.strip()
            ls.append(line)
        passscore = eval(candidatels2[len(ls)-1][1])
        candils = []
        for i in candidatels2:
            if eval(i[1])>=0.7*passscore:
                candils.append(i[0])
            else:
                break
        with open(dislodgeout,'w') as out:
            for i in ls:
                if i in candils:
                    out.write(i+'\n')
          

## main
argv = sys.argv[1]
# prepare hmm profile
with open(argv,'r') as infile:
    gene = []
    for line in infile:
        if line.startswith('#'):
            continue
        line = line.strip()
        if line == '':
            continue
        line = line.split()
        fasta = line[0]
        genename = line[1]
        if genename not in gene:
            gene.append(genename)
        else:
            print(" Warning ".center(40,"-"))
            print('[Warning]:There is repetitive gene name '+ genename+', please confirm they are same, otherwise only the first will be effective')
            print('Please enter "c" for continue or "e" for exit:')
            i = input()
            if i == 'c':
                continue
            elif i == 'e':
                exit()
            else:
                print('Unrecognized input, this script will exit.')
                exit()
        ids = line[2]
        id_list = ids.split(',')
        Athfasta = sourcepath+'Arabidopsis_thaliana.fa'
        genenamefa = temppath+genename+'.fa'
        extract_seq(id_list,Athfasta,genenamefa)
        print('extract sequence: OK')
        alignin = genenamefa
        alignout = temppath+genename+'.a.fa' 
        align(alignin,alignout)
        print('align: OK')
        fa = alignout
        hmmfile = temppath+genename+'.hmm'
        hmmpre(fa,hmmfile,genename)
        print('Hmm profile: OK')
# hmmscan , blastp trace, dislodge
with open(argv,'r') as infile:
    for line in infile:
        if line.startswith('#'):
            continue
        testhomopass = False
        line = line.strip()
        if line == '':
            continue
        line = line.split()
        fasta = line[0]
        genename = line[1]
        ids = line[2]
        id_list = ids.split(',')
        if testhomo(id_list,genename) == False:
            print(" Warning ".center(40,"-"))
            print("[Warning]:This family ("+genename+") has not passed the homology test, please cheak if the members are of the same family or complete.")
            break
        hmmfile = temppath+genename+'.hmm'
        fas = sourcepath+fasta
        hmmout = temppath+fasta+'_'+genename+'.can.out'
        hmmtblout = temppath+fasta+'_'+genename+'.can.tblout'
        hmm(hmmfile,hmmout,hmmtblout,fas)
        candidatels,candidatels2 = analyse(hmmtblout)
        canfa = temppath+genename+fasta+'.can.fa'
        extract_seq(candidatels,fas,canfa)
        query = canfa
        blastout = temppath+fasta+'_'+genename+'.can.blast'
        targetnum = len(id_list)
        blast(query,referencedb,blastout,targetnum)
        traceout = temppath+fasta+'_'+genename+'.trace.txt'
        trace(blastout,id_list,traceout)
        dislodgeout = outpath+fasta+'_'+genename+'.txt'
        dislodge(traceout,candidatels2,dislodgeout)
        finalfa = outpath+fasta+'_'+genename+'.fa'
        with open(dislodgeout,'r') as refile:
            ls2 = []
            for line in refile:
                line = line.strip()
                ls2.append(line)
            extract_seq(ls2,fas,finalfa)


end_time = datetime.datetime.now()
print('')
print(' Well Done '.center(40,'='))
print(str(end_time-start_time).center(40))
