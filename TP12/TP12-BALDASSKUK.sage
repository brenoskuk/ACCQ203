print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP12 : CRYPTANALYSE ET CRYPTOGRAPHIE A BASE DE RESEAUX                      #
# *************************************************************************** #
# *************************************************************************** #
""")

# CONSIGNES
#
# Les seules lignes a modifier sont annoncee par "Code pour l'exercice"
# indique en commmentaire et son signalees
# Ne changez pas le nom des variables
#
# CONSEILS
#
# Ce modele vous sert a restituer votre travail. Il est deconseille d'ecrire
# une longue suite d'instruction et de debugger ensuite. Il vaut mieux tester
# le code que vous produisez ligne apres ligne, afficher les resultats et
# controler que les objets que vous definissez sont bien ceux que vous attendez.
#
# Vous devez verifier votre code en le testant, y compris par des exemples que
# vous aurez fabrique vous-meme.
#


reset()
print("""\
# ****************************************************************************
# SAC A DOS (NE PAS TRAITER LA QUESTION 2)
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

b = [356,278,417,27,132,464,521]
s = 1287

# Code pour l'EXERCICE

M = matrix(ZZ,   [[1,0,0,0,0,0,0,0],
                 [0,1,0,0,0,0,0,0],
                 [0,0,1,0,0,0,0,0],
                 [0,0,0,1,0,0,0,0],
                 [0,0,0,0,1,0,0,0],
                 [0,0,0,0,0,1,0,0],
                 [0,0,0,0,0,0,1,0],
                 [b[0],b[1],b[2],b[3],b[4],b[5],b[6],-s]])

x = M.T.LLL()[0]

# # Affichage des resultats

print "Le message est"

print x

sp = 0
for i in range(len(b)):
    sp += x[i]*b[i]  

assert sp == s

reset()
print("""\
# ****************************************************************************
# ATTAQUE DE WIENER
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

N1 = 65946239999
e1 = 22022476093
N2 = 65946239999
e2 = 10865199773

# Code pour l'EXERCICE

# N1 = N2...

N = N1

X = int(N^(1/4))
Y = int(N^(3/4))

b1_1 = [e1*X, Y]
b1_2 = [e2*X, Y]
b2 = [N*X, 0]

A1 = matrix(ZZ, [b1_1, b2]).T
A2 = matrix(ZZ, [b1_2, b2]).T


# we expect that the smallest element returned by LLL is of help
b_1_lll = A1.T.LLL()[0]
b_2_lll = A2.T.LLL()[0]


s1 = (b_1_lll[1]/Y)
k1 = -((b_1_lll[0] - s1*e1*X)/(N*X))


s2 = -(b_2_lll[1]/Y)
k2 = -((b_2_lll[0] - s2*e2*X)/(N*X))


sigma1 = -int((b_1_lll[0])/(k1*X)-1)
sigma2 = -int((b_2_lll[0])/(k2*X)-1)


Pol.<x> = PolynomialRing(RR)

pol1 = x^2 + sigma1*x + N1
pol2 = x^2 + sigma2*x + N1


p1,q1 = pol1.roots()
p1 = int(p1[0])
q1 = int(q1[0])
p2,q2 = pol2.roots()
p2 = int(p2[0])
q2 = int(q2[0])


# # Affichage des resultats

print "Using weiner attack we find:\nCandidates for p,q given e1 = ", p1,',',q1
print "Candidates for p,q given e2 = ", p2,',',q2
print "Since when p and q are similar we may extract them using the attack, we prioritize e1."


reset()
print("""\
# ****************************************************************************
# METHODE DE COPPERSMITH
# ****************************************************************************
""")

from math import e as euler
import numpy as np
# Donnees de l'enonce de l'exercice

Pol.<x> = PolynomialRing(ZZ)

# Code pour l'EXERCICE


def Coppersmith(f,N,B=-1):

    P = f.parent()
    x = P.gen()
    d = f.degree()

    m = ceil(log(N)/d)
    
    if B == -1:
        B = ceil(N^(1/d)/(2*euler))
    
    coeffs = f.coefficients(sparse = False)
    # f = x^d + ad-1*x^d-1 + ... + a0

    M = np.zeros((2*d,2*d))
    for i in range(d):
        M[i,i] = N^m*B^i

    for i in range(d):
        for j in range(d):
            M[i+j,d+j] = coeffs[i]*B^(i+j)


    for i in range(d):
        M[d+i][d+i] = B^(d+i)

    M = Matrix(ZZ,M)

    V = M.T.LLL()

    v = list(V[0])

    h = sum([v[i]*x^i/(B^i) for i in range(len(v))])

    R = h.roots(QQ,multiplicities=False)
    rep = [r for r in R if (abs(r)<B and mod(f(r),N)==0)]
    return rep

# # Affichage des resultats

p=(x+1)*(x-2)*(x-3)*(x-29)
print Coppersmith(p,10000,30)



reset()
import numpy as np
from math import e as euler
print("""\
# ****************************************************************************
# MESSAGES STEREOTYPES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

bin=BinaryStrings()
N = 42564360034887861127
Pol.<x> = PolynomialRing(ZZ)
PolmodN.<y> = PolynomialRing(Integers(N))
e = 3
c = 12843085802751039909

# Code pour l'EXERCICE

def Coppersmith(f,N,B=-1):

    P = f.parent()
    x = P.gen()
    d = f.degree()

    m = ceil(log(N)/d)
    
    if B == -1:
        B = ceil(N^(1/d)/(2*euler))
    
    coeffs = f.coefficients(sparse = False)
    # f = x^d + ad-1*x^d-1 + ... + a0

    M = np.zeros((2*d,2*d))
    for i in range(d):
        M[i,i] = N^m*B^i

    for i in range(d):
        for j in range(d):
            M[i+j,d+j] = coeffs[i]*B^(i+j)


    for i in range(d):
        M[d+i][d+i] = B^(d+i)

    M = Matrix(ZZ,M)

    V = M.T.LLL()

    v = list(V[0])

    h = sum([v[i]*x^i/(B^i) for i in range(len(v))])

    R = h.roots(QQ,multiplicities=False)
    rep = [r for r in R if (abs(r)<B and mod(f(r),N)==0)]
    return rep



text = "jj/mm:XX"
# we generate te approximate message
mt = mod( ZZ(str(bin.encoding(text)), base=2),N)

ZN = IntegerModRing(N)
# we know that or approximate message differs from c in a few places
Pol.<x> = PolynomialRing(ZN)
f = (mt + x)^e - c

# Coopersmith fails due to N being too big... found no fix
#Coppersmith(f,N,B=-1)
roots = f.roots(multiplicities=False)
r = roots[0]
mr = mt + r


# to binary
mrb=ZZ(mr).str(2)

m = '0'*(8*ceil(len(mrb)/8)-len(mrb))+mrb
decoded_m = bin(m).decoding()


# # Affichage des resultats

print "Ce jour la, le message est"
print decoded_m

reset()
print("""\
# ****************************************************************************
# ALGORITHME DE BABAI (NE TRAITER QUE LA QUESTION 1)
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE

def Babai(B,t):
    b = vector(t)
    # calculate orthogonal basis
    G , _ = B.T.gram_schmidt()
    
    G = G.T
    
    # get dimensions
    n = G.ncols()
    
    for i in range(n-1,-1,-1):
        bt = G.column(i)
        u = b.dot_product(bt)/bt.dot_product(bt)
        b = b - round(u)*B.column(i)
    return t-b

# # Affichage des resultats

B = matrix(ZZ, [
                [1,2,0,0,5],
                [1,6,1,0,10],
                [3,0,0,1,-2],
                [1,0,2,0,-7],
                [1,0,2,7,3]
               ])

Bl = B.LLL()

v1 = vector([1,1,1,1,3])

v2 = vector([1,1,1,1,4])

vt1 = Babai(Bl,v1)

vt2 = Babai(Bl,v2)

print 'For reduced base \n', Bl
print 'for vector v1 ', v1
print 'closest vector to v1 in reduced base is ', vt1
print 'for vector v2 ', v2
print 'closest vector to v2 in reduced base is ', vt2


reset()
print("""\
# ****************************************************************************
# CRYPTOSYSTEME GGH
# ****************************************************************************
""")


def StringToAscii(text):
    return [ZZ(ord(s)) for s in list(text)]
def AsciiToString(listascii):
    return "".join([str(unichr(a)) for a in listascii])
def StringToInts(text):
    bin = BinaryStrings()
    return [ZZ(str(t),base=2) for t in list(bin.encoding(text))]
def IntsToString(listints):
    bin = BinaryStrings()
    return bin("".join([str(t) for t in listints])).decoding()

# Donnees de l'enonce de l'exercice

EE= matrix(ZZ,40, {(32, 32): 1, (22, 13): -2, (23, 8): -3,
       (15, 31): -1, (22, 37): -1,
       (13, 5): -4, (38, 20): 2, (4, 12): 3, (19, 22): -2, (15, 5): 2, (11,
       32): -1, (11, 10): 3, (1, 11): -4, (12, 33): 1, (0, 15): 1, (33, 17): 1,
       (7, 19): -1, (11, 1): -2, (7, 27): 3, (19, 32): -4, (22, 10): 2, (31,
       39): -4, (34, 9): 2, (36, 17): 2, (18, 17): 1, (14, 6): -2, (23, 14): 3,
       (23, 34): 2, (12, 11): -3, (0, 21): -3, (27, 22): -2, (4, 29): -3, (23,
       5): 1, (4, 6): -2, (24, 7): 2, (5, 38): -2, (33, 13): -1, (9, 35): 3,
       (18, 36): 1, (22, 5): 1, (24, 25): 3, (34, 31): 2, (6, 34): -3, (23,
       33): -4, (20, 37): -1, (38, 12): 2, (33, 0): -1, (4, 32): 3})
AA=10*identity_matrix(40)+EE
HH=AA.hermite_form()
cc = vector([-2, 0, 2, 0, 0, 1, -1, -1, -3, 0, 0, 2, -1, 13, 7, 2, 0, 2, 27, 2, 1,
       17, -2, 899, 50, 15, 11, 1098, 7, 2, -1, 10, -1, 2, 156, 15, 42, 8,
       525748584, 37])


# Code pour l'EXERCICE

def Babai(B,t):
    b = vector(t)
    # calculate orthogonal basis
    G , _ = B.T.gram_schmidt()
    
    G = G.T
    
    # get dimensions
    n = G.ncols()
    
    for i in range(n-1,-1,-1):
        bt = G.column(i)
        u = b.dot_product(bt)/bt.dot_product(bt)
        b = b - round(u)*B.column(i)
    return t-b


# solution
mt = Babai(AA, cc)
UU = AA.T.inverse()*HH.T
m = (UU.inverse().T*mt)

# # Affichage des resultats

print "given function didn't work... in any case the message to be decoded is \n", m


