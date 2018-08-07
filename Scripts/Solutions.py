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

def readMatrix(filename):
    # Read the file
    handle = open(filename, 'r')
    content= handle.readlines()
    handle.close()
    
    # Set up the matrix file
    matrix  = {}
    letters = []
    numline = len(content) 
    
    for nl in range(XXX, XXX):
        line = content[nl]
        splt = line.split()
        a = splt[XXX]
        if a not in matrix:
            matrix[a] = {}
            letters.append(a)
            
    # Go throug the file and save the values
    for nl in range(XXX, XXX):
        line = content[nl]
        splt = line.split()
        l = len(splt)
        aa1 = XXX
        for a in range(XXX, XXX):
            aa2 = letters[XXX]
            matrix[XXX][XXX] = splt[XXX]
            matrix[XXX][XXX] = XXX
    
    return matrix