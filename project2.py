def sDNA_to_mRNA(sDNA):
    mRNA = ''
    for c in sDNA:
        if c == 'T':
            mRNA += 'U'
        else:
            mRNA += c
    return mRNA


def rev_comp(RNA0):
    RNA1 = ''
    for c in RNA0[::-1]:
        if c == 'A':
            RNA1 += 'U'
        if c == 'U':
            RNA1 += 'A'   
        if c == 'G':
            RNA1 += 'C'   
        if c == 'C':
            RNA1 += 'G'
    return RNA1


def Ui_Tei(gsiRNA, psiRNA):
    count = 0
    d = 0
    e = 0
    f = True
    if gsiRNA[0] in 'AU':
        count+=1
    if psiRNA[0] in 'GC':
        count+=1
    if (psiRNA[:7].count('A') + psiRNA[:7].count('U')) >= 5:
        count+=1
    for i in range(21):
        if gsiRNA[i] in 'GC':
            d+=1
        else:
            d=0
        if psiRNA[i] in 'GC':
            e+=1
        else:
            e=0        
        if d>9 or e>9:
            f = False
            break
    if f == True:
        count+=1
    return count


def Reynolds(gsiRNA, psiRNA):
    global if_loop
    if_loop = False
    count = 0
    if 30<=((psiRNA.count('G') + psiRNA.count('C') + gsiRNA.count('G') + gsiRNA.count('C'))/42*100)<=52:
        count+=1
    if (psiRNA[16:].count('A') + psiRNA[16:].count('U'))>=3:
        count+=1
    for i in range(8):
        term = rev_comp(gsiRNA[i:i+6])
        if term in gsiRNA[i+8:]:
            if_loop = True
            break
    for i in range(8):
        term = rev_comp(psiRNA[i:i+6])
        if term in psiRNA[i+8:]:
            if_loop = True
            break    
    if if_loop == False:
        count+=1
    if psiRNA[18]=='A':
        count+=1
    if psiRNA[2] == 'A':
        count+=1    
    if psiRNA[9]=='U':
        count+=1  
    if psiRNA[18] in 'AU':
        count+=1    
    if psiRNA[12]!='G':
        count+=1  
    return count


def Amar(gsiRNA, psiRNA):
    count = 0
    if (gsiRNA[14:21].count('A') + gsiRNA[14:21].count('U'))<(gsiRNA[1:8].count('A') + gsiRNA[1:8].count('U')):
        count+=1
    if psiRNA[0] in 'GC':
        count+=1
    if psiRNA[0] != 'U':
        count+=1 
    if psiRNA[5] == 'A':
        count+=1  
    if psiRNA[18] in 'AU':
        count+=1 
    if psiRNA[18] != 'G':
        count+=1     
    return count

maxres = 0
DNA = ''
with open (input('Path to .fasta file on your computer with DNA: '))  as file_in:
    for line in file_in:
        if line.startswith('>'):
            continue
        else:
            DNA += line
            if DNA.endswith('\n'):
                DNA = DNA[:len(DNA)-1]
DNA = DNA.upper()
if (DNA.count('A')+ DNA.count('T')+DNA.count('G')+DNA.count('C'))==len(DNA):
    mRNA = sDNA_to_mRNA(DNA)
else:
    mRNA = ''
    print('Incorrect input')
l_a = []
l_b = []
l_mRNA_site = []
l_siRNA_guide = []
l_siRNA_passenger = []
for i in range(2, len(mRNA)-20):
    gsiRNA = rev_comp(mRNA[i-2:i+19])
    psiRNA = mRNA[i:i+21]
    R = Reynolds(gsiRNA, psiRNA)
    if if_loop == False:
        res = Ui_Tei(gsiRNA, psiRNA) + R + Amar(gsiRNA, psiRNA)
    else:
        res = 0
    if res > maxres:
        l_a = []
        l_b = []
        l_mRNA_site = []
        l_siRNA_guide = []
        l_siRNA_passenger = []
        maxres = res
        a = i-1
        b = i+19
        mRNA_site = mRNA[i-2:i+19]
        siRNA_guide = gsiRNA
        siRNA_passenger = psiRNA
        efficacy = res/18*100
        l_a.append(a)
        l_b.append(b)
        l_mRNA_site.append(mRNA_site)
        l_siRNA_guide.append(siRNA_guide)
        l_siRNA_passenger.append(siRNA_passenger)             
    elif res == maxres:
        maxres = res
        a = i-1
        b = i+19
        mRNA_site = mRNA[i-2:i+19]
        siRNA_guide = gsiRNA
        siRNA_passenger = psiRNA
        efficacy = res/18*100
        l_a.append(a)
        l_b.append(b)
        l_mRNA_site.append(mRNA_site)
        l_siRNA_guide.append(siRNA_guide)
        l_siRNA_passenger.append(siRNA_passenger)    
with open (input('Path to .txt file on your computer for output: '), 'w') as file_out:
    for i in range(len(l_a)):
        file_out.write(f"mRNA_site [{l_a[i]}, {l_b[i]}] (5'->3'):")       
        file_out.write('\n')
        file_out.write(l_mRNA_site[i])
        file_out.write('\n')
        file_out.write("siRNA_guide (5'->3'):")
        file_out.write('\n')
        file_out.write(l_siRNA_guide[i])
        file_out.write('\n')
        file_out.write("siRNA_passenger (5'->3'):")
        file_out.write('\n')
        file_out.write(l_siRNA_passenger[i])
        file_out.write('\n')
        file_out.write(f'efficacy: {efficacy}% ({maxres}/18)')
        file_out.write('\n')
        file_out.write('\n')
        file_out.write('\n')
    file_out.write(str(len(l_a)))