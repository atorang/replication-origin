def PatternCount(Text, Pattern):
    "count pattern in the text"
    count=0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)]==Pattern:
            count+=1
    return count
####################################################################################################
def FrequentWords(Text, k):
    'Find the most frequent k-mer'
    FrequentPatterns=[]
    Count=[]
    for i in range(len(Text)-k+1):
        Pattern=Text[i:k+i]
        Count.append(PatternCount(Text,Pattern))
    maxCount=max(Count)
    for i in range(len(Text)-k+1):
        if Count[i]==maxCount:
            FrequentPatterns.append(Text[i:k+i])
    FrequentPatterns=set(FrequentPatterns)  #remove duplicates
    return FrequentPatterns
########################################################################################
def ReverseComplement(Text):
    "Generate Reverse Complement of the string"
    RC=''
    Complement={'A':'T','C':'G','G':'C','T':'A'}
    Text=reversed(Text)
    for i in Text:
        RC=RC+Complement[i]
    return RC
#######################################################################################
def PatternMatching(Pattern, Genome):
    'A collection of integers specifying all starting positions where Pattern appears as a substring of Genome'
    Position = []
    for i in range(len(Genome) - len(Pattern)+1):
        if Genome[i:i + len(Pattern)] == Pattern:
            Position.append(i)
    return Position
#######################################################################################
def SymbolToNumber(symbol):
    S={'A':0,'C':1,'G':2,'T':3}
    return S[symbol]
#######################################################################################
def NumberToSymbol(index):
    S = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in S.keys():
        if S[i]==index:
            return i
#######################################################################################
def PatternToNumber(Pattern):
    if len(Pattern)==0:
        return 0
    symbol=Pattern[-1]
    Prefix=Pattern[:-1]
    return 4*PatternToNumber(Prefix)+SymbolToNumber(symbol)
########################################################################################
def NumberToPattern(index, k):
    if k == 1:
        return NumberToSymbol(index)
    Quotient=divmod(index, 4)
    prefixIndex=Quotient[0]
    r=Quotient[1]
    symbol=NumberToSymbol(r)
    PrefixPattern=NumberToPattern(prefixIndex,k-1)
    return PrefixPattern+symbol
#########################################################################################
def ComputingFrequencies(Text, k):
    FrequencyArray=[]
    for i in range(4**k):
        FrequencyArray.append(0)
    for i in range(len(Text)-k+1):
        Pattern=Text[i:i+k]
        j=PatternToNumber(Pattern)
        FrequencyArray[j]=FrequencyArray[j]+1
    return FrequencyArray
#########################################################################################
def FindingFrequentWordsBySorting(Text , k):
    FrequentPatterns=[]
    Index=[]
    Count=[]
    for i in range(len(Text)-k+1):
        Pattern=Text[i:i+k]
        Index.append(PatternToNumber(Pattern))
        Count.append(1)
    SortedIndex=sorted(Index)
    for i in range(len(Text)-k+1):
        if SortedIndex[i]==SortedIndex[i-1]:
            Count[i] = Count[i-1]+1
    maxCount=max(Count)
    for i in range(len(Text)-k+1):
        if Count[i] == maxCount:
            Pattern=NumberToPattern(SortedIndex[i], k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
#########################################################################################
def ClumpFinding(Genome,k,L,t):
    #We defined a k-mer as a "clump" if it appears many times within a short interval of the genome.
    # More formally, given integers L and t, a k-mer Pattern forms an (L, t)-clump inside a (longer) string Genome
    #  if there is an interval of Genome of length L in which this k-mer appears at least t times.
    FrequentPatterns=[]
    Clump=[]
    for i in range(4**k):
        Clump.append(0)
    for i in range(len(Genome)-L+1):
        Text=Genome[i:i+L]
        FrequencyArray=ComputingFrequencies(Text, k)
        for index in range(4**k):
            if FrequencyArray[index]>= t:
                Clump[index]=1
    for i in range(4**k):
        if Clump[i]== 1:
            Pattern=NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
##########################################################################################
def BetterClumpFinding(Genome, k, L, t):
    FrequentPatterns=[]
    Clump=[]
    for i in range(4**k):
        Clump.append(0)
    Text=Genome[0:L]
    FrequencyArray=ComputingFrequencies(Text, k)
    for i in range(4**k):
        if FrequencyArray[i] >= t:
            Clump[i]=1
    for i in range(1,len(Genome)-L+1):
        FirstPattern=Genome[i-1:i-1+k]
        index=PatternToNumber(FirstPattern)
        FrequencyArray[index]=FrequencyArray[index]-1
        LastPattern=Genome[i + L-k: i+L]
        index =PatternToNumber(LastPattern)
        FrequencyArray[index]=FrequencyArray[index] + 1
        if FrequencyArray[index] >= t:
            Clump[index]=1
    for i in range(4**k):
        if Clump[i] == 1:
            Pattern=NumberToPattern(i,k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
####################################################################################
def Skew(Genome):
    set={'C':-1,'G':1,'A':0,'T':0}
    skew=[0]
    for i in range(len(Genome)):
        skew.append(skew[i]+set[Genome[i]])
    return skew
#############################################################################
def MinimumSkew(Genome):
    #Find a position in a genome where the skew diagram attains a minimum
    skew=Skew(Genome)
    Mins=[]
    mins=min(skew)
    for i in range(len(skew)):
        if skew[i]==mins:
            Mins.append(i)
    return Mins
###################################################################
def HammingDistance(p,q):
    #The number of mismatches between strings p and q is called the Hamming distance between these strings
    count = 0
    for i in range(len(p)):
        if p[i]!=q[i]:
            count+=1
    return count
###################################################################
def ApproximatePatternMatching(Genome,Pattern,d):
    'A collection of integers specifying the position of all approximate occurrences of a pattern in a string'
    Position = []
    for i in range(len(Genome) - len(Pattern)+1):
        if HammingDistance(Genome[i:i + len(Pattern)],Pattern)<=d:
            Position.append(i)
    return Position
####################################################################
def ApproximatePatternCount(Text, Pattern,d):
    "count pattern in the text"
    count=0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i + len(Pattern)],Pattern)<=d:
            count+=1
    return count
####################################################################
def ImmediateNeighbors(Pattern):
    Neighborhood=[Pattern]
    for i in range(len(Pattern)):
        symbol=Pattern[i]
        for j in ['A','C','G','T']:
            if symbol!= j:
                Neighbor=Pattern[0:i]+j+Pattern[i+1:]
                Neighborhood.append(Neighbor)
    return Neighborhood
#####################################################################
def Neighbors(Pattern, d):
    # generate all k-mers of Hamming distance less than d from Pattern.
    if d == 0:
        return Pattern
    if len(Pattern)== 1:
        return ['A', 'C', 'G', 'T']
    Neighborhood=[]
    SuffixNeighbors=Neighbors(Pattern[1:], d)
    for i in range(len(SuffixNeighbors)):
        if HammingDistance(Pattern[1:], SuffixNeighbors[i]) < d:
            for j in ['A','C','G','T']:
                Neighborhood.append(j+SuffixNeighbors[i])
        else:
            Neighborhood.append(Pattern[0] + SuffixNeighbors[i])
    return Neighborhood
#####################################################################
def IterativeNeighbors(Pattern, d):
    # generate all k-mers of Hamming distance exactly d from Pattern.
    Neighborhood=[Pattern]
    for j in range(d):
        for i in range(len(Neighborhood)):
            Neighborhood.extend(ImmediateNeighbors(Neighborhood[i]))
            Neighborhood=set(Neighborhood)
            Neighborhood=list(Neighborhood)
    return Neighborhood
#####################################################################
def FrequentWordsWithMismatches(Text, k, d):
    FrequentPatterns=[]
    Neighborhoods=[]
    Index=[]
    Count=[]
    for i in range(len(Text)-k+1):
        Neighborhoods.extend(Neighbors(Text[i:i+k], d))
    NeighborhoodArray=Neighborhoods
    for i in range(len(Neighborhoods)):
        Pattern=NeighborhoodArray[i]
        Index.append(PatternToNumber(Pattern))
        Count.append(1)
    SortedIndex=sorted(Index)
    for i in range(len(Neighborhoods)-1):
        if SortedIndex[i] == SortedIndex[i + 1]:
            Count[i + 1]=Count[i] + 1
    maxCount=max(Count)
    for i in range(len(Neighborhoods)):
        if Count[i] == maxCount:
            Pattern=NumberToPattern(SortedIndex[i], k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
#############################################################################
def FrequentWordsWithMismatchesAndReverseComplements (Text, k, d):
    FrequentPatterns=[]
    Neighborhoods=[]
    Index=[]
    Count=[]
    Reverse=[]
    for i in range(len(Text)-k+1):
        Neighborhoods.extend(Neighbors(Text[i:i+k], d))
    for i in range(len(Neighborhoods)):
        Reverse.append(ReverseComplement(Neighborhoods[i]))
    NeighborhoodArray=Neighborhoods
    NeighborhoodArray.extend(Reverse)
    for i in range(len(Neighborhoods)):
        Pattern=NeighborhoodArray[i]
        Index.append(PatternToNumber(Pattern))
        Count.append(1)
    SortedIndex=sorted(Index)
    for i in range(len(Neighborhoods)-1):
        if SortedIndex[i] == SortedIndex[i + 1]:
            Count[i + 1]=Count[i] + 1
    maxCount=max(Count)
    print('Max count=',maxCount)
    for i in range(len(Neighborhoods)):
        if Count[i] == maxCount:
            Pattern=NumberToPattern(SortedIndex[i], k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
############################################################
def MotifEnumeration(Dna, k, d):
    Patterns=[]
    k_mer=[]
    for i in range(len(Dna[0])-k+1):
        k_mer.append(Dna[0][i:i+k])
    for i in range(len(k_mer)):
        Neighborhoods=Neighbors(k_mer[i],d)
        if Neighborhoods==k_mer[i]:
            for n in range(len(Dna)):
                if ApproximatePatternCount(Dna[n], Neighborhoods,d)==0:
                    break
            else:
                Patterns.append(Neighborhoods)
        else:
            for j in range(len(Neighborhoods)):
                for n in range(len(Dna)):
                    if ApproximatePatternCount(Dna[n], Neighborhoods[j],d)==0:
                        break
                else:
                    Patterns.append(Neighborhoods[j])
    Patterns=set(Patterns)
    return list(Patterns)
##############################################################
def Count(Motifs):
    count=[]
    if len(Motifs) != 1:
        for i in range(len(Motifs[0])):
            A = 0;C = 0;G = 0; T = 0
            for j in range(len(Motifs)):
                if Motifs[j][i] == 'A': A += 1
                if Motifs[j][i] == 'C': C += 1
                if Motifs[j][i] == 'G': G += 1
                if Motifs[j][i] == 'T': T += 1
            count.append([A, C, G, T])
        count=[list(x) for x in zip(*count)]
    else:
        for j in range(len(Motifs[0])):
            A = 0;C = 0;G = 0; T = 0
            if Motifs[0][j] == 'A': A += 1
            if Motifs[0][j] == 'C': C += 1
            if Motifs[0][j] == 'G': G += 1
            if Motifs[0][j] == 'T': T += 1
            count.append([A, C, G, T])
        count = [list(x) for x in zip(*count)]
    return count
##############################################################
def Profile(Motifs):
    count=Count(Motifs)
    profile=count
    for i in range(4):
        for j in range(len(Motifs[0])):
            profile[i][j]=count[i][j]/len(Motifs)
    return profile
##############################################################
def Consensus(Motifs):
    set = {0:'A', 1:'C', 2:'G', 3:'T'}
    consensus=''
    if len(Motifs[0])!=1:
        for i in range(len(Motifs[0])):
            A = 0; C = 0; G = 0; T = 0
            for j in range(len(Motifs)):
                if Motifs[j][i]=='A': A+=1
                if Motifs[j][i]=='C': C+=1
                if Motifs[j][i]=='G': G+=1
                if Motifs[j][i]=='T': T+=1
            S=[A,C,G,T]
            Max=max(S)
            for n in range(4):
                if S[n]==Max:
                    consensus=consensus+set[n]
                    break
    else: print('Motifs has one string')
    return consensus
##############################################################
def Score(Motifs):
    consensus=Consensus(Motifs)
    score=0
    for i in range(len(Motifs)):
        score+=HammingDistance(consensus,Motifs[i])
    return score
##############################################################
def D_Pattern_Motifs(Pattern,Motifs):
    D=0
    for i in range(len(Motifs)):
        D+=HammingDistance(Pattern,Motifs[i])
    return D
###############################################################
def Motifs_Pattern_Dna(Pattern, Dna):
    motifs=[]
    for i in range(len(Dna)):
        distance = []
        for j in range(len(Dna[0])-len(Pattern)+1):
            distance.append(HammingDistance(Dna[i][j:j+len(Pattern)],Pattern))
        Min=min(distance)
        for j in range(len(Dna[0]) - len(Pattern) + 1):
            if Min==distance[j]:
                motifs.append(Dna[i][j:j+len(Pattern)])
                break
    return motifs
###########################################################
def D_Pattern_Text(Pattern,Text):
    k=len(Pattern)
    D = []
    for j in range(len(Text)-k+1):
        D.append(HammingDistance(Pattern,Text[j:j+k]))
    distance=min(D)
    return distance
###########################################################
def Motif(Pattern, Text):
    k = len(Pattern)
    distance=D_Pattern_Text(Pattern, Text)
    motif=''
    for j in range(len(Text) - k + 1):
        if HammingDistance(Pattern, Text[j:j + k])==distance:
            motif=Text[j:j + k]
            break
    return motif
###########################################################
def D_Pattern_Dna(Pattern,Dna):
    distance=0
    for i in range(len(Dna)):
        distance+=D_Pattern_Text(Pattern,Dna[i])
    return distance
###############################################################
def MedianString(Dna, k):
    # find a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern,
    # the same task that the Equivalent Motif Finding Problem is trying to achieve.
    # We call such a k-mer a median string for Dna.
    distance=k*len(Dna)
    Patterns=[]
    Median=''
    for index in range(4**k):
        Patterns.append(NumberToPattern(index,k))
    for i in range(len(Patterns)):
        if distance > D_Pattern_Dna(Patterns[i], Dna):
            distance=D_Pattern_Dna(Patterns[i], Dna)
            Median=Patterns[i]
    return Median
###############################################################
def Probability(Pattern, Profile):
    set = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    probability=1
    k=len(Pattern)
    for j in range(k):
        probability = probability * Profile[set[Pattern[j]]][j]
    return probability
###############################################################
def ProfileMostProbableKmer(Text, k, profile):
    prob=[]
    Kmer=[]
    for i in range(len(Text)-k+1):
        pattern=Text[i:i+k]
        probability=Probability(pattern,profile)
        prob.append(probability)
    maxProb=max(prob)
    for i in range(len(Text)-k+1):
        if prob[i]==maxProb:
            Kmer=Text[i:i+k]
            break
    return Kmer
################################################################
def GreedyMotifSearch(Dna, k, t):
    BestMotifs=[]
    for i in range(len(Dna)):
        BestMotifs.append(Dna[i][:k])
    for i in range(len(Dna[0])-k+1):
        Motif1=Dna[0][i:i+k]
        Motifs=[]
        Motifs.append(Motif1)
        for j in range(1,t):
            Motifi=ProfileMostProbableKmer(Dna[j],k,Profile(Motifs))
            Motifs.append(Motifi)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs=Motifs
    return BestMotifs
##############################################################
def Profile_Pseudocounts(Motifs,pseudocount):
    count=Count(Motifs)
    profile=count
    for i in range(4):
        for j in range(len(Motifs[0])):
            profile[i][j]=(count[i][j]+pseudocount)/(len(Motifs)+pseudocount*4)
    return profile
##################################################################
def GreedyMotifSearchWithPseudocounts(Dna, k, t, pseudocount):
    BestMotifs = []
    for i in range(len(Dna)):
        BestMotifs.append(Dna[i][:k])
    for i in range(len(Dna[0]) - k + 1):
        Motif1 = Dna[0][i:i + k]
        Motifs = []
        Motifs.append(Motif1)
        for j in range(1, t):
            Motifi = ProfileMostProbableKmer(Dna[j], k, Profile_Pseudocounts(Motifs,pseudocount))
            Motifs.append(Motifi)
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
##################################################################
def DistanceBetweenPatternAndStrings(Pattern, Dna):
    distance=0
    for i in range(len(Dna)):
        distance+=D_Pattern_Text(Pattern,Dna[i])
    return distance
##################################################################
def Motifs_Profile_Dna(Profile,Dna):
    motifs=[]
    k=len(Profile[0])
    for i in range(len(Dna)):
        prob=[]
        for j in range(len(Dna[0])-k+1):
            Pattern=Dna[i][j:j+k]
            prob.append(Probability(Pattern,Profile))
        Max=max(prob)
        for j in range(len(Dna[0])):
            if prob[j]==Max:
                motifs.append(Dna[i][j:j+k])
                break
    return motifs
##################################################################
import random
def RandomizedMotifSearch(Dna, k, t):
    motifs=[]
    for i in range(t):
        x=random.randint(0,len(Dna[0])-k)
        motifs.append(Dna[i][x:x+k])
    BestMotifs=motifs
    Profile=[]
    while True:
        Profile=Profile_Pseudocounts(motifs,1)
        motifs=Motifs_Profile_Dna(Profile, Dna)
        if Score(motifs) < Score(BestMotifs):
            BestMotifs=motifs
        else:
            return BestMotifs
##################################################################
def RandomizedMotifSearch_nTime(Dna, k, t, n):
    score = []
    motifs = []
    for i in range(n):
        m = RandomizedMotifSearch(Dna, k, t)
        score.append(Score(m))
        motifs.append(m)
    Min = min(score)
    Finalmotifs = []
    for i in range(n):
        if score[i] == Min:
            Finalmotifs = motifs[i]
    return Finalmotifs
##################################################################
def ProfileRandomlyGeneratedKmer(Text,profile,k):
    prob=[]
    for i in range(len(Text)-k+1):
        pattern=Text[i:i+k]
        prob.append(Probability(pattern,profile))
    Sum=sum(prob)
    for i in range(len(Text)-k+1):
        prob[i]=prob[i]/Sum
    Range = [0]
    for i in range(len(Text) - k+1):
        Range.append(prob[i]+Range[i])
    Range.append(1)
    Random=random.uniform(0,1)
    RandomPattern = ''
    for i in range(len(Text) - k+1):
        if Random>=Range[i] and Random<=Range[i+1]:
            RandomPattern=Text[i:i+k]
    return RandomPattern
##################################################################
def GibbsSampler(Dna, k, t, N):
    motifs = []
    for i in range(t):
        x = random.randint(0, len(Dna[0]) - k)
        motifs.append(Dna[i][x:x + k])
    BestMotifs = motifs.copy()
    for j in range(N):
        i=random.randint(0,t-1)
        motifs_i=motifs[0:i]+motifs[i+1:]
        Profile = Profile_Pseudocounts(motifs_i, 1)
        motifs[i]=ProfileRandomlyGeneratedKmer(Dna[i],Profile,k)
        if Score(motifs) < Score(BestMotifs):
            BestMotifs=motifs.copy()
    return BestMotifs
###################################################################
def GibbsSampler_nTime(Dna, k, t, N, repeat):
    score = []
    motifs = []
    for i in range(repeat):
        m = GibbsSampler(Dna,k,t,N)
        score.append(Score(m))
        motifs.append(m)
    Min = min(score)
    Finalmotifs = []
    for i in range(repeat):
        if score[i] == Min:
            Finalmotifs = motifs[i]
    return Finalmotifs