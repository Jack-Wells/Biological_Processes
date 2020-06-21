#Packages
import numpy as np
import matplotlib.pyplot as plt
import random
from sklearn.preprocessing import normalize
from Bio import Entrez
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
import Bio
Entrez.email = ''

#Part 1
#Implimentation of Needleman-Wunsch algorithm with backtracking
def bNW(s1, s2, match, diff, d_diff, d_none):
    #This part produces the Needleman-Wunsch matrix
    NW = np.zeros((len(s1)+1,len(s2)+1))
    
    for i in range(0,len(s1)+1):
        NW[i][0] = diff * i
    
    # Fill out first row
    for j in range(0,len(s2)+1):
        NW[0][j] = diff * j

    for i in range(1,len(s1)+1):
        for j in range(1,len(s2)+1):

            if s1[i-1] == s2[j-1]:
                Match = NW[i-1][j-1] + match
            else:
                Match = NW[i-1][j-1] + d_diff 
                
            left = NW[i-1][j] + diff
            top = NW[i][j-1] + diff
            # Record the maximum score from the three possible scores calculated above
            NW[i][j] = max(Match, left, top)
    #Backtracking to find the al       
    x,y = NW.shape
    
    i = x-1
    j = y-1

    align1 = ''
    align2 = ''
    
    Path = []
    Path.append([i,j])
    while i > 0 and j > 0: # end touching the top or the left edge
        score_current = NW[i][j]
        score_diagonal = NW[i-1][j-1]
        score_up = NW[i][j-1]
        score_left = NW[i-1][j]
        
        # Check to figure out which cell the current score was calculated from
        if s1[i-1] == s2[j-1]:
            d = match
        else:
            d = d_diff
        if score_current == score_up + diff:
            align1 += '-'
            align2 += s2[j-1]
            j -= 1
            Path.append([i,j])
        elif score_current == score_left + diff:
            align1 += s1[i-1]
            align2 += '-'
            i -= 1
            Path.append([i,j])
        elif score_current == score_diagonal + d:
            align1 += s1[i-1]
            align2 += s2[j-1]
            i -= 1
            j -= 1
            Path.append([i,j])

    #To go to top left cell
    while i > 0:
        align1 += s1[i-1]
        align2 += '-'
        i -= 1
        Path.append([i,j])
    while j > 0:
        align1 += '-'
        align2 += s2[j-1]
        j -= 1
        Path.append([i,j])

    #As we traversed the score matrix from the bottom right, our two sequences will be reversed.
    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2, NW, Path

sequence1 = 'CAAAGACCTGAAGAGCCAGTGGACTCCACCCCACTTTCTGGTCTGACCAATT'
sequence2 = 'ACCACACTCTCTGGGCTGACCAATTACAGCGCTTCTACAGAACTGAAGACTCC'

s1,s2,m,p = bNW(sequence1,sequence2,1,-0.5,0,None)
print(s1)
print(s2)
print(m)
print(p)

score = 0 
for i in range(len(s1)):
    if s1[i] == s2[i]:
        score += 1
    elif s1[i] == '-' and s2[i] == '-':
        score += 0
    elif s1[i] == '-' and s2[i] != '-':
        score += -0.5
    elif s2[i] == '-' and s1[i] != '-':
        score += -0.5
print(score)

#Image of the path taken by backtracking
image = np.ones((len(sequence1)+1,len(sequence2)+1))
for i in p:
    image[i[0],i[1]] = 255
plt.imshow(image)
plt.show()

#Task 2
#Acessing the gene sequences
handle = Entrez.efetch(db="nucleotide", id="NM_001030.6",rettype="fasta", retmode="xml")
record = Entrez.read(handle,"fasta")
handle.close()
record=record[0] #record is an array of dictionaries, get the first entry
human = record['TSeq_sequence']

handle = Entrez.efetch(db="nucleotide", id="NM_027015.4",rettype="fasta", retmode="xml")
record = Entrez.read(handle,"fasta")
handle.close()
record=record[0] #record is an array of dictionaries, get the first entry
mouse = record['TSeq_sequence']

handle = Entrez.efetch(db="nucleotide", id="NM_001251997.1",rettype="fasta", retmode="xml")
record = Entrez.read(handle,"fasta")
handle.close()
record=record[0] #record is an array of dictionaries, get the first entry
chimp = record['TSeq_sequence']

handle = Entrez.efetch(db="nucleotide", id="XM_026928592.1",rettype="fasta", retmode="xml")
record = Entrez.read(handle,"fasta")
handle.close()
record=record[0] #record is an array of dictionaries, get the first entry
catfish = record['TSeq_sequence']

handle = Entrez.efetch(db="nucleotide", id="XM_006030956.2",rettype="fasta", retmode="xml")
record = Entrez.read(handle,"fasta")
handle.close()
record=record[0] #record is an array of dictionaries, get the first entry
alligator = record['TSeq_sequence']

animals = [human, mouse, chimp, catfish, alligator]


M = np.ones((len(animals),len(animals)))

for i in range(len(animals)):
    for j in range(len(animals)):
        if i == j:
            M[i][j] = 0
        else:
            a1,a2,a,b = bNW(animals[i],animals[j],1,-1,-1, -np.inf)
            S = 0 #matches
            for k in range(len(a1)):
                if a1[k] == a2[k]:
                    S += 1
                elif a1[k] == '-' and a2[k] == '-':
                    S += np.inf
                elif a1[k] != a2[k]:
                    S -= 1
            N = len(a1) #backtrack length
            d = (N-S)/N #distance
            M[i][j] = float(d)
print(M)

#Builds a lower triangular matrix 
low_tri = [] 
for i in range(len(animals)):
    low_tri.append(list(M[i,0:i+1]))
    
for i in range(len(low_tri)):
    for j in range(i+1):
        low_tri[i][j] = float(low_tri[i][j])

#produces a distance matrix for the tree
distance = DistanceMatrix(['human', 'mouse', 'chimp', 'catfish', 'alligator'],low_tri)
DTC = DistanceTreeConstructor()
tree = DTC.nj(distance)
#Gives a description and lengths of the branches of the tree
print(tree)
#draws the tree
Phylo.draw(tree)

#part 4 - Toggle switch

lamda = 0.25
alpha = 250
ki = 200
ka = 400
hi = 2
ha = 20

def dg_e(g_a, g_b):
    value = alpha * 1/(1 + (g_b/ki)**hi) - lamda*g_a
    return value

def dg_l(g_c, g_a, g_d):
    value = alpha * 1/(1 + (ka/g_a)**ha) * 1/(1 + (g_d/ki)**hi) - lamda*g_c
    return value
	
dt = 1
sigma = 10

#Inintail conditions
g_a = []
g_b = []
g_c = []
g_d = []
g_e = []
g_f = []
    
current_g_a = 300
current_g_b = 300
current_g_c = 0
current_g_d = 0
current_g_e = 0
current_g_f = 0
    
g_a.append(current_g_a)
g_b.append(current_g_b)
g_c.append(current_g_c)
g_d.append(current_g_d)
g_e.append(current_g_e)
g_f.append(current_g_f)

for t in range(100):
    new_g_a = current_g_a + dt*dg_e(current_g_a,current_g_b) + np.sqrt(sigma)*np.random.normal()
    new_g_b = current_g_b + dt*dg_e(current_g_b,current_g_a) + np.sqrt(sigma)*np.random.normal()
    
    new_g_c = current_g_c + dt*dg_l(current_g_c,current_g_a,current_g_d) + np.sqrt(sigma)*np.random.normal()
    new_g_d = current_g_d + dt*dg_l(current_g_d,current_g_a,current_g_c) + np.sqrt(sigma)*np.random.normal()
    new_g_e = current_g_e + dt*dg_l(current_g_e,current_g_b,current_g_f) + np.sqrt(sigma)*np.random.normal()
    new_g_f = current_g_f + dt*dg_l(current_g_f,current_g_b,current_g_e) + np.sqrt(sigma)*np.random.normal()

    g_a.append(new_g_a)
    g_b.append(new_g_b)
    g_c.append(new_g_c)
    g_d.append(new_g_d)
    g_e.append(new_g_e)
    g_f.append(new_g_f)

    current_g_a = new_g_a
    current_g_b = new_g_b
    current_g_c = new_g_c
    current_g_d = new_g_d
    current_g_e = new_g_e
    current_g_f = new_g_f

#Plots of the gene activation
t = range(0,101)
fig, ax = plt.subplots(4, figsize=(8,8))
fig.suptitle('Activation of downstream genes')
plt.subplot(2,2,1)
plt.title('gC')
plt.xlabel('Time')
plt.ylabel('Activation')
plt.plot(t,g_c,color='red')

plt.subplot(2,2,2)
plt.title('gD')
plt.xlabel('Time')
plt.ylabel('Activation')
plt.plot(t,g_d,color='blue')

plt.subplot(2,2,3)
plt.title('gE')
plt.xlabel('Time')
plt.ylabel('Activation')
plt.plot(t,g_e,color='green')

plt.subplot(2,2,4)
plt.xlabel('Time')
plt.ylabel('Activation')
plt.title('gF')
plt.plot(t,g_f,color='black')

plt.subplots_adjust(wspace=0.4, hspace=0.4)
plt.show()

dt = 1
sigma = 10
data = np.zeros((400,101,6))

for i in range(400):
    #Inintail conditions
    g_a = []
    g_b = []
    g_c = []
    g_d = []
    g_e = []
    g_f = []
    
    current_g_a = 300
    current_g_b = 300
    current_g_c = 0
    current_g_d = 0
    current_g_e = 0
    current_g_f = 0
    
    g_a.append(current_g_a)
    g_b.append(current_g_b)
    g_c.append(current_g_c)
    g_d.append(current_g_d)
    g_e.append(current_g_e)
    g_f.append(current_g_f)

    for t in range(100):
        new_g_a = current_g_a + dt*dg_e(current_g_a,current_g_b) + np.sqrt(sigma)*np.random.normal()
        new_g_b = current_g_b + dt*dg_e(current_g_b,current_g_a) + np.sqrt(sigma)*np.random.normal()
    
        new_g_c = current_g_c + dt*dg_l(current_g_c,current_g_a,current_g_d) + np.sqrt(sigma)*np.random.normal()
        new_g_d = current_g_d + dt*dg_l(current_g_d,current_g_a,current_g_c) + np.sqrt(sigma)*np.random.normal()
        new_g_e = current_g_e + dt*dg_l(current_g_e,current_g_b,current_g_f) + np.sqrt(sigma)*np.random.normal()
        new_g_f = current_g_f + dt*dg_l(current_g_f,current_g_b,current_g_e) + np.sqrt(sigma)*np.random.normal()

        g_a.append(new_g_a)
        g_b.append(new_g_b)
        g_c.append(new_g_c)
        g_d.append(new_g_d)
        g_e.append(new_g_e)
        g_f.append(new_g_f)

        data[i][t+1] = [new_g_a, new_g_b, new_g_c, new_g_d, new_g_e, new_g_f]

        current_g_a = new_g_a
        current_g_b = new_g_b
        current_g_c = new_g_c
        current_g_d = new_g_d
        current_g_e = new_g_e
        current_g_f = new_g_f


#Producing samples of the data
samples = np.zeros((400,6))

for i in range(400):
    time = int(random.random()*101) #random time to sample the data
    samples[i] = data[i][time]

print(samples)


#Emmbedding the 6 dimensional data in 3D
laplacian = np.zeros((400,400)) #laplacian matrix

for i in range(400):
    for j in range(400):
        if i == j:
            laplacian[i,j] = 0
        elif i != j:
            #Formula from lectures
            laplacian[i,j] = np.exp(-(1/2) * (np.linalg.norm(samples[i] - samples[j])/sigma)**2)

norm = normalize(laplacian) #Normalizing the laplacian

value, vector = np.linalg.eig(norm) #Calculates the eigenvalues and eigenvectors
value = list(value)
vector = list(vector)
largest = []

#Takes the largest 3 eigenvalues and their corresponding eigenvectors
for i in range(3):
    large = value.index(max(value))
    largest.append(list(vector[large]))
    del value[large]
    del vector[large]

print(largest) 

