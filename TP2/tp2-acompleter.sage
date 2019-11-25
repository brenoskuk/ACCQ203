print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP2 : ALGEBRE LINEAIRE SUR UN ANNEAU PRINCIPAL                              #
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
# MISE SOUS FORME NORMALE D'HERMITE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
        [-2,  3,  3,  1],
        [ 2, -1,  1, -3],
        [-4,  0, -1, -4]])

A1 = random_matrix(ZZ, 7, 8, algorithm='echelonizable', rank=3)

# Code pour l'EXERCICE

def MyHNF(A):
    """
    Forme normale d'Hermite selon votre code
    """
    m = A.nrows()
    n = A.ncols()
    H = copy(A)
    U = identity_matrix(n)
    # ECRIVEZ VOTRE CODE ICI
    assert(H-A*U==0)
    return H,U

def SageHNF(A):
    """
    Forme normale d'Hermite de SAGE avec la convention francaise :
    Les vecteurs sont ecrits en colonne.
    """
    m = A.nrows()
    n = A.ncols()
    Mm = identity_matrix(ZZ,m)[::-1,:]
    Mn = identity_matrix(ZZ,n)[::-1,:]
    AA = (Mm*A).transpose()
    HH, UU = AA.hermite_form(transformation=True)
    H = (HH*Mm).transpose()*Mn
    U = UU.transpose()*Mn
    assert(H-A*U==0)
    return H,U

H,  U  = MyHNF(A)
HH, UU = SageHNF(A)

test = False # Test a ecrire

# # Affichage des resultats

print "\n$ Question 4"
print "La matrice A = "
print A
print "a pour forme normale d'Hermite H="
print H
print "et matrice de transformation U="
print U
print "\n$ Question 5"
print "D'apres SageMath, la matrice A a pour forme normale d'Hermite H="
print HH
print "et matrice de transformation U="
print UU
print "\n$ Question 6"
print "Les deux fonctions coincident-elles sur une centaine d'exemples ?"
print test

reset()
print("""\
# ****************************************************************************
# MISE SOUS FORME NORMALE DE SMITH
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

X1 = matrix(ZZ,2,3,[
        [4, 7, 2],
        [2, 4, 6]])

X2 = matrix(ZZ,3,3,[
        [-397, 423, 352],
        [   2,  -3,   1],
        [-146, 156, 128],
])

PolQ.<xQ> = PolynomialRing(QQ)
AQ = matrix(PolQ,3,[
            [xQ + 1,  2,     -6],
            [     1, xQ,     -3],
            [     1,  1, xQ - 4]])
Pol2.<x2> = PolynomialRing(FiniteField(2))
AF2 = matrix(Pol2,3,[
            [x2 + 1,  2,     -6],
            [     1, x2,     -3],
            [     1,  1, x2 - 4]])
Pol3.<x3> = PolynomialRing(FiniteField(3))
AF3 = matrix(Pol3,3,[
            [x3 + 1,  2,     -6],
            [     1, x3,     -3],
            [     1,  1, x3 - 4]])
Pol5.<x5> = PolynomialRing(FiniteField(5))
AF5 = matrix(Pol5,3,[
            [x5 + 1,  2,     -6],
            [     1, x5,     -3],
            [     1,  1, x5 - 4]])

# Code pour l'EXERCICE

def MySNF(A):
    """
    Forme normale de Smith selon votre code
    """
    m = A.nrows()
    n = A.ncols()
    H = copy(A)
    L = identity_matrix(m)
    U = identity_matrix(n)
    # ECRIVEZ VOTRE CODE ICI
    assert(H-L*A*U==0)
    return H,L,U

H1, L1, U1 = MySNF(X1)
H2, L2, U2 = MySNF(X2)

HQ, _, _ = MySNF(AQ)
HF2, _, _ = MySNF(AF2)
HF3, _, _ = MySNF(AF3)
HF5, _, _ = MySNF(AF5)

test = False # Test a ecrire

# # Affichage des resultats

print "\n$ Question 4"
print "La matrice X1 = "
print X1
print "a pour forme normale de Smith H1="
print H1
print "et matrice de transformation L1="
print L1
print "et matrice de transformation U1="
print U1
print "La matrice X2 = "
print X2
print "a pour forme normale de Smith H2="
print H2
print "et matrice de transformation L2="
print L2
print "et matrice de transformation U2="
print U2

print "\n$ Question 5"
print "La forme normale de Smith sur Q est "
print AQ
print "La forme normale de Smith sur F2 est "
print AF2
print "La forme normale de Smith sur F3 est "
print AF3
print "La forme normale de Smith sur F5 est "
print AF5

print "\n$ Question 6"
print "Votre fonction coincide avec celle de Sage ?"
print test


reset()
print("""\
# ****************************************************************************
# RESOLUTION DE SYSTEME LINEAIRE HOMOGENE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

X = matrix(ZZ,[
      [ -2,  3,  3,  1],
      [  2, -1,  1, -3],
      [ -4,  0, -1, -4]])

# Code pour l'EXERCICE

L =[] # liste des vecteurs d'une base

# # Affichage des resultats

print "Le systeme a pour racine le module engendre par"
print L

reset()
print("""\
# ****************************************************************************
# IMAGE D'UNE MATRICE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A  = matrix(ZZ, [
           [ 15,  8, -9, 23,  -9],
           [ 22, 22,  7, -8,  20],
           [ 21, 18, -1, -7,  -3],
           [  3, -1,  0, 12, -16]])


# Code pour l'EXERCICE

test = false # test a ecrire

# # Affichage des resultats

print "L'image de"
print A
print "est-elle egale a ZZ^4 ?"
print test



reset()
print("""\
# ****************************************************************************
# RESOLUTION DE SYSTEME LINEAIRE NON-HOMOGENE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

X1 = matrix(ZZ,[
           [ -6,  12,  6],
           [ 12, -16, -2],
           [ 12, -16, -2]])

b1 = vector(ZZ,[ -6, 4, 4])

PolF5.<x> = PolynomialRing(GF(5))

X2 = matrix(PolF5,[
           [ x + 1, 2,     4],
           [     1, x,     2],
           [     1, 1, x + 1]])

b2 = vector(PolF5,[ 3*x+2, 0, -1])

# Code pour l'EXERCICE

z1 = vector(ZZ,3)
H1 = []

z2 = vector(PolF5,3)
H2 = []

z3 = vector(ZZ,3)
H3 = []

# # Affichage des resultats

print "Une solution particuliere de X1*z1 = b1 est"
print z1
print "les solutions du systeme homogene sont engendres par"
print H1
print "Une solution particuliere de X2*z2 = b2 est"
print z2
print "les solutions du systeme homogene sont engendrees par"
print H2
print "Une solution particuliere du systeme 3 est"
print z3
print "les solutions du systeme homogene sont engendres par"
print H3



reset()
print("""\
# ****************************************************************************
# STRUCTURE DU QUOTIENT
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
              [ -630,   735,   0,   735, -630],
              [ 1275, -1485, -15, -1470, 1275],
              [  630,  -630,   0,  -630,  630]])

# Code pour l'EXERCICE

reponse = "A ecrire"

# # Affichage des resultats

print "La structure de Z^3/N est"
print reponse

reset()
print("""\
# ****************************************************************************
# FACTEURS INVARIANTS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
              [ -630,   735,   0,   735, -630],
              [ 1275, -1485, -15, -1470, 1275],
              [  630,  -630,   0,  -630,  630]])


# Code pour l'EXERCICE

rang = -1 
fact_inv = []
reponse = "A ecrire"

# # Affichage des resultats

print "Le rang de Z^3 / N est"
print rang
print "Les facteurs invariants sont"
print fact_inv
print "Exposants ?"
print reponse

