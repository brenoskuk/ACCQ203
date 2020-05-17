print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP10 : INVARIANTS DE SIMILITUDE ET LFSR                                     #
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
import numpy as np
print("""\
# ****************************************************************************
# INVARIANTS DE SIMILITUDE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

PolQQ.<x>=PolynomialRing(QQ)

M1 =  matrix(QQ,[[-1841, -10363, 22304, 108021, -243809 ],
[1366, 7695, -16535, -80130, 180869 ],[
-1072, -6088, 13069, 63408, -143144 ],[
506, 2951, -6298, -30700, 69343 ],[
82, 502, -1061, -5214, 11788]])

M2 = matrix(QQ,[[570, 1652, -8251, 3807, 34007 ],[
-178, -522, 2666, -1196, -10988 ],[
540, 1573, -7866, 3622, 32430 ],[
-42, -118, 580, -275, -2387 ],[
135, 393, -1967, 905, 8109]])

M3 = matrix(QQ,[[64, -300, 924, -228, 3168 ],[
-80, 404, -1232, 304, -4224 ],[
35, -175, 543, -133, 1848 ],[
-15, 75, -231, 61, -792 ],[
-20, 100, -308, 76, -1052]])


# Code pour l'EXERCICE

# Q1

# Mp = M - xI

Mp1 = M1 - x*matrix(QQ,np.eye(5))
Mp2 = M2 - x*matrix(QQ,np.eye(5))
Mp3 = M3 - x*matrix(QQ,np.eye(5))

Delta1 =  Mp1.smith_form()[0]
Delta2 =  Mp2.smith_form()[0]
Delta3 =  Mp3.smith_form()[0]

Invariants1 = []
Invariants2 = []
Invariants3 = []


Polcar1 = 1
Polcar2 = 1
Polcar3 = 1
for i in range(5):
    if Delta1[i][i] != 1:
        Invariants1.append(Delta1[i][i])
    if Delta2[i][i] != 1:
        Invariants2.append(Delta2[i][i])
    if Delta3[i][i] != 1:
        Invariants3.append(Delta3[i][i])
    Polcar1 = Delta1[i][i]*Polcar1
    Polcar2 = Delta2[i][i]*Polcar2
    Polcar3 = Delta3[i][i]*Polcar3

Polmin1 = Invariants1[len(Invariants1) - 1]
Polmin2 = Invariants2[len(Invariants2) - 1]
Polmin3 = Invariants3[len(Invariants3) - 1]


# # Affichage des resultats

print "La forme normale de M1 est\n", Delta1
print "La forme normale de M2 est\n", Delta2
print "La forme normale de M3 est\n", Delta3

print "\nLes invariants de similitude de M1 sont\n", Invariants1
print M1.rational_form(format='invariants')
print "Les invariants de similitude de M2 sont\n", Invariants2
print M2.rational_form(format='invariants')
print "Les invariants de similitude de M3 sont\n", Invariants3
print M3.rational_form(format='invariants')

print "\nLe polynôme minimal de M1 est ",Polmin1
print "et d'après SageMath", M1.minimal_polynomial()
print "Le polynôme minimal de M2 est ",Polmin2
print "et d'après SageMath", M2.minimal_polynomial()
print "Le polynôme minimal de M3 est ",Polmin3
print "et d'après SageMath", M3.minimal_polynomial()

print "\nLe polynôme caractéristique de M1 est ",Polcar1
print "et d'après SageMath", M1.characteristic_polynomial()
print "\nLe polynôme caractéristique de M2 est ",Polcar2
print "et d'après SageMath", M2.characteristic_polynomial()
print "\nLe polynôme caractéristique de M3 est ",Polcar3
print "et d'après SageMath", M3.characteristic_polynomial()




print("""\
# ****************************************************************************
# FORME DE FROBENUIS D'UNE MATRICE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

# idem exercice précédent



# Code pour l'EXERCICE

# creates the frobenius matrix from a list of invariants
def get_Frobenius(P):
    c = 0
    
    F = np.zeros((5,5))
    
    for i in range(len(P)):
        p_list = P[i].list()
        count = 0
        r = len(p_list) - 1
        for j in range(r):
            if(j + 1 < r):
                F[c+j + 1][c+j] = 1
            F[c+j][c+r-1] = p_list[j]
            count += 1
            
        c = c + count
    return F
    
    
Frob1 = get_Frobenius(Invariants1)
Frob2 = get_Frobenius(Invariants2)
Frob3 = get_Frobenius(Invariants3)

# # Affichage des resultats

print "La forme de Frobenius de M1 est\n", Frob1
print "vérification Sagemath \n", M1.rational_form()
print "La forme de Frobenius de M2 est\n", Frob2
print "vérification Sagemath \n", M2.rational_form()
print "La forme de Frobenius de M3 est\n", Frob3
print "vérification Sagemath \n", M3.rational_form()


print("""\
# ****************************************************************************
# FORME DE JORDAN D'UNE MATRICE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

# idem exercice précédent

# Code pour l'EXERCICE

# given list of invariants, get irreductibles and create Jordan form
def get_Jordan(P):
    c = 0
    
    J = np.zeros((5,5))
    
    # decompose polynomial
    Pp = []
    for poly in P:
        factors = poly.factor()
        for fac in factors:
            Pp.append(fac)
    
    # sort 
    Pp = sorted(Pp)

    for i in range(len(Pp)):
        deg = Pp[len(Pp) - i - 1][1]
        fac = Pp[len(Pp) - i - 1][0]    
        #print fac, ' ** ', deg
        count = 0
        for j in range(deg):
            if(j - 1 >= 0):
                J[c+j - 1][c+j] = 1
            J[c+j][c+j] = -1*fac[0]
            count += 1
        c = c + count         
    return J


Jordan1 = get_Jordan(Invariants1)
Jordan2 = get_Jordan(Invariants2)
Jordan3 = get_Jordan(Invariants3)

# # Affichage des resultats

print "La forme de Jordan de M1 est\n", Jordan1
print "vérification Sagemath \n", M1.jordan_form()
print "La forme de Jordan de M2 est\n", Jordan2
print "vérification Sagemath \n", M2.jordan_form()
print "La forme de Jordan de M3 est\n", Jordan3
print "vérification Sagemath \n", M3.jordan_form()


reset()
print("""\
# ****************************************************************************
# BERLEKAMP-MASSEY
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

F2 = GF(2)
suite = [F2(1),0,1,1,0,0,0,0,1]
suite = suite+suite+suite+suite
ans = berlekamp_massey(suite)


# Code pour l'EXERCICE


def myBerlekampMassey(r):
    Fq = r[0].parent()
    PolFq.<x> = PolynomialRing(Fq)
    c = PolFq(1)
    cp = PolFq(1)
    l = 0
    dp = PolFq(1)
    m = -1
    d = r[0]
    sigma = 0
    for k in range(len(r)):
        # discrepancy
        
        # create list of coefficients of c 
        coeffs = [0]*(len(r))
        clist = PolFq(c).list()
                    
        for j in range(len(clist)):
            coeffs[j] = clist[j]
        
        # sum 
        sigma = Fq(0)
        for i in range(1,l+1):
            #print ('l=',l,'\nc = ',c ,'\ncoeffs = ', coeffs, '\ni - 1', i - 1, '\nk-i = ',k-i)
            
            sigma += Fq(coeffs[i])*r[k-i]
        
        d = r[k] + sigma
        if d != 0:
            t = c
            c = c - d*(dp^(-1))*cp*x^(k-m)
            if l <= k/2:
                l = k + 1 - l
                cp = t
                dp = d
                m = k    
    return c,l
    
    

# # Affichage des resultats

q=2
Fq = FiniteField(q)
for _ in range(10):
    t=randint(1,10)
    r = [Fq.random_element() for _ in range(t) ]
    r = r+r+r+r
    p1 = berlekamp_massey(r)
    p2,_ = myBerlekampMassey(r)
    
    print 'p1 = ', p1
    print 'p2 = ', p2
    print '\n'
    if p1.reverse()- p2 !=0:
        print "Erreur"
reset()
print("""\
# ****************************************************************************
# ORBITES D'UN LFRS
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice
# Code pour l'EXERCICE
print "NE PAS TRAITER"
# # Affichage des resultats

reset()
print("""\
# ****************************************************************************
# PERIODES D'UN LFRS
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

F2 = GF(2)
PolF2.<x> = PolynomialRing(F2)
chi1 = 1 + x + x^2 + x^3 + x^4
chi2 = 1 + x + x^2 + x^4
chi3 = 1 + x^3 + x^6

F5 = GF(5)
PolF5.<y> = PolynomialRing(F5)
chi4 = 3 + y^2
chi5 = 3 + 3*y + y^2


# Code pour l'EXERCICE

def BSGS(chi):
    PolFq = chi.parent()
    x = PolFq.gen()
    Fq = PolFq.base_ring()
    n = len(chi.list()) - 1

    m = ceil(sqrt(n))
    T = []
    e = Fq(1)
    for j in range(m):
        T.append(e)
        e = e*chi
        if j > 0 and PolFq(e) == Fq(1):
            return j
    print 'm = ', m, ' T = ', T, '\n'
        
    h = chi^(-m)
    gamma = h
    for i in range(1,m):
    
        for j in range(0,m):
            if gamma == T[j]:
                return i*m + j
        gamma = gamma*h
            

# # Affichage des resultats

print "Période de chi1",  BSGS(chi1)
print "Période de chi2",  BSGS(chi2)
print "Période de chi3",  BSGS(chi3)
print "Période de chi4",  BSGS(chi4)
print "Période de chi5",  BSGS(chi5)

reset()
print("""\
# ****************************************************************************
# SERIES GENERATRICES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

# Code pour l'EXERCICE


print "NE PAS TRAITER"

# # Affichage des resultats
