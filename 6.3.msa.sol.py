#!/usr/bin/env python

import sys
import re
class Seq:
    def __init__(self):
        self.id="x"
    def __init__(self):
        self.seq="y"
        self.features="z"

class Node:
    def __init__(self):
        self.name=""
        self.distance=float (-1)
        self.bootstrap=float (-1)
        self.left=0
        self.right=0
        self.parent=0

def declare_new_tree_node(nodes,nn):
    nodes[nn]=Node()
    return (nn, nn+1)

def readmatrix(filename):
  handle = open(filename, "r")
  content = handle.readlines()
  handle.close()

  matrix = {}
  letters = []
  nlines=len (content)


  for ln in range (0,nlines-1):
      line=content[ln]
      splitted = line.split()
      a=splitted[0]
      if a not in matrix:
          matrix[a] = {}
          letters.append(a)

  for ln in range (0,nlines-1):
      line=content[ln]
      splitted = line.split()
      l=len(splitted)
      aa1=splitted[0]
      for a in range (1,l):
          aa2=letters[a-1]
          matrix[aa1][aa2]=splitted[a]
          matrix[aa2][aa1]=splitted[a]

  return matrix

def align (seq,Igroup, Jgroup, matrix, gep):

    lenI=len(seq[Igroup[0]])
    lenJ=len(seq[Jgroup[0]])

    smat = [[0 for x in range(lenJ+1)] for y in range(lenI+1)]
    tb   = [[0 for x in range(lenJ+1)] for y in range(lenI+1)]

    for i in range (0, lenI+1):
      smat[i][0]=i*gep
      tb[i][0]=1

    for j in range (0, lenJ+1):
      smat[0][j]=j*gep
      tb[0][j]=-1


    for i in range (1, lenI+1):
      for j in range (1, lenJ+1):
        s=float(0)
        nsub=float(0)
        for ni in range (0,len(Igroup)):
            for nj in range (0, len(Jgroup)):

               a1=seq[Igroup[ni]][i-1]
               a2=seq[Jgroup[nj]][j-1]
               if a1!='-' and a2!='-':
                  s+=int(matrix[a1.upper()][a2.upper()])
                  nsub+=1
        if (nsub>0):
           s/=nsub

        Sub=smat[i-1][j-1]+s
        Del=smat[i][j-1]+gep
        Ins=smat[i-1][j]+gep

        if Sub>Del and Sub >Ins:
            smat[i][j]=Sub
            tb  [i][j]=0
        elif Del>Ins:
            smat[i][j]=Del
            tb[i][j]=-1
        else:
            smat[i][j]=Ins
            tb[i][j]=1


    #print "Optimal Score: %d\n"%(int(smat[lenI][lenJ]))
    i=lenI
    j=lenJ
    lenA=0
    alnI=[]
    alnJ=[]
    new_seq={}
    for ni in range (0,len (Igroup)):
        new_seq[Igroup[ni]]=""
    for nj in range (0,len (Jgroup)):
        new_seq[Jgroup[nj]]=""

    while ((i==0 and j==0)!=1):
        if (tb[i][j]==0):
            i-=1
            j-=1
            for ni in range (0,len (Igroup)):
                new_seq[Igroup[ni]]+=seq[Igroup[ni]][i]
            for nj in range (0,len (Jgroup)):
                new_seq[Jgroup[nj]]+=seq[Jgroup[nj]][j]

        elif (tb[i][j]==-1):
            j-=1
            for ni in range (0,len (Igroup)):
                new_seq[Igroup[ni]]+="-"
            for nj in range (0,len (Jgroup)):
                new_seq[Jgroup[nj]]+=seq[Jgroup[nj]][j]

        elif (tb[i][j]==1):
            i-=1
            for ni in range (0,len (Igroup)):
                new_seq[Igroup[ni]]+=seq[Igroup[ni]][i]
            for nj in range (0,len (Jgroup)):
                new_seq[Jgroup[nj]]+="-"

        lenA+=1

    for ni in range (0,len (Igroup)):
        seq[Igroup[ni]]=new_seq[Igroup[ni]][::-1]
    for nj in range (0,len (Jgroup)):
        seq[Jgroup[nj]]=new_seq[Jgroup[nj]][::-1]
    return seq

def node2splits (N,nodes,seq,matrix,gep):
    list=[]
    if nodes[N].name!="":
        list.append(nodes[N].name)
    else:
        left_list=[]
        rigt_list=[]
        if nodes[N].left:
            left_list=node2splits(nodes[N].left, nodes,seq,matrix,gep)
        if nodes[N].right:
            right_list=node2splits(nodes[N].right, nodes,seq,matrix,gep)

        if (len (right_list)>0 and len(left_list)>0):
                seq=align (seq, right_list, left_list,matrix,gep)

        list=left_list+right_list

    return list

def readseq (fname):
    seqlist={}
    record=[]
    nrec=-1
    inseq=0
    with open(fname) as f:
      for line in f:
        if re.match ( r'^>', line):
            nrec+=1
            record.append(Seq())
            mobj=re.match ( r'^>(\S*)\s*(.*)', line)

            if (mobj):
                record[nrec].id=mobj.group(1)
                record[nrec].features=mobj.group(2)
            inseq=0
        else :
            if inseq==0 :
                inseq=1
                record[nrec].seq=line
            else:
                cstring=record[nrec].seq+line
                record[nrec].seq=cstring

    seqlist={}
    for x in range (0,nrec+1):
        record[x].seq=re.sub (r'[ \n\t\r-]',"",record[x].seq)
        seqlist[record[x].id]=record[x].seq

    return seqlist



tree=""

with open (sys.argv[1]) as f:
    for line in f:
        tree+=line

tree=re.sub (r'[ \n\t\r]',"",tree)

#get the tree
nn=0
nodes={}
(nodes, nn)=newick2nodes(tree)

seq=readseq(sys.argv[2])
matrix=readmatrix (sys.argv[3])

node2splits(1,nodes, seq, matrix, -4)
for s in seq.keys():
    sys.stdout.write (">%s\n%s\n"%(s, seq[s]))
