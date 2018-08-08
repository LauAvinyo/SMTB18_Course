def readFastaFile(file):
    sequences = {}
    
    with open(file) as f:
        for line in f:
            # If this line is name line
            if re.match(r'^>', line):
                matchObject = re.match(r'^>(\S*)\s*(.*)', line)
                
                if (matchObject):
                    name = matchObject.group(1)
                    sequences[name] = ''
                inseq = 0
            # Otherwise
            else:
                if inseq == 0:
                    inseq = 1
                    sequences[name] = line.strip()
                else:
                    sequences[name] += line.strip()
    return sequences
    

def generateMatrix(alignment):
    totResidues = 0
    totSubtitut = 0
    
    # Initialize to zero the total residue counting
    freqr = {}
    for res in alp:
        freqr[res] = 0
    
    # Initialize to zero the subtition matrix
    matrix = {}
    for res1 in alp:
        matrix[res1] = {}
        for res2 in alp:
            matrix[res1][res2] = 0
            
    
    # Count how many time each residue is found
    # Mind the Gap!
    for s in alignment:
        for aa in alignment[s]:
            if aa!='-': 
                freqr[aa] += 1
                totResidues += 1
    
    # Get the frequency of each residue
    for aa in freqr:
        freqr[aa] /= totResidues
    
    # Since dictionaries are not good to numerate in a fix manner
    # Lets transform the values into a tuple
    
    sequences = tuple(alignment[i] for i in alignment)
    nseq = len(sequences)
    
    for i in range(0, nseq):
        for j in range(i+1, nseq):
            seqI = sequences[i]
            seqJ = sequences[j]
            if len(seqI) != len(seqJ): return "ERROR!!"
            
            for aa in range(0, len(seqI)):
                aaI = seqI[aa]
                aaJ = seqJ[aa]
                if aaI != '-' and aaJ != '-':
                    matrix[aaI][aaJ] += 1
                    matrix[aaJ][aaI] += 1
                    totSubtitut += 2
    
    # mathematical transformation 
    for aa1 in alp:
        for aa2 in alp:
            num = matrix[aa1][aa2] / totSubtitut
            den = freqr[aa1]*freqr[aa2]
            matrix[aa1][aa2] = math.log10(num/den) * 10
        
    return matrix


seqI = 'THEFASTCAT'
seqJ = 'THEFATCAT'

def pairAlignment(seqI, seqJ, matrix):
    
    # Check that there are no gaps in this sequences
    # I am paranoid
    seqI.replace('-', '')
    seqJ.replace('-', '')

    # Get the length of the sequences
    lenI=len(seqI)
    lenJ=len(seqJ)

    # Penalty for putting a gap
    gep=-4

    
    # Initiallize to zero the matrices
    smat = [[0 for x in range(lenJ+1)] for y in range(lenI+1)]
    tb   = [[0 for x in range(lenJ+1)] for y in range(lenI+1)]

    # Base cases
    for i in range (0, lenI+1):
        smat[i][0]=0
        tb[i][0]=-2

    for j in range (0, lenJ+1):
        smat[0][j]=0
        tb[0][j]=-2

    # Fill the table
    bscore=0
    for i in range (1, lenI+1):
        for j in range (1, lenJ+1):
            if seqI[i-1]!='-' and seqJ[j-1]!='-':
                s=int(matrix[seqI[i-1]][seqJ[j-1]])
            else:
                s=0

            Sub=smat[i-1][j-1]+s
            Del=smat[i][j-1]+gep
            Ins=smat[i-1][j]+gep

            if Sub>Del and Sub >Ins and Sub >0:
                smat[i][j]= Sub
                tb  [i][j]= 0
            elif Del>Ins and Del >0:
                smat[i][j]= Del
                tb[i][j]=-1
            else:
                smat[i][j]= Ins
                tb[i][j]=1

            if smat[i][j]>bscore:
                bscore=smat[i][j]

    # Traceback
    alnI=''
    alnJ=''
    while (tb[i][j]!=-2):
        if (tb[i][j]==0):
            i-=1
            j-=1
            alnI += seqI[i]
            alnJ += seqJ[j]
        elif (tb[i][j]==-1):
            j-=1
            alnI += '-'
            alnJ += seqJ[j]
        elif (tb[i][j]==1):
            i-=1
            alnI += seqI[i]
            alnJ += '-'


    # This is used to reverse the sequences
    alnI=alnI[::-1]
    alnJ=alnJ[::-1]

    return alnI, alnJ, bscore