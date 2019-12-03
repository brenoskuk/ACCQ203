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
    Forme normale d'Hermite avec Matrice de passage U
    """
    m = A.nrows()
    n = A.ncols()
    H = copy(A)
    U = identity_matrix(n)
    # CODE:
    l = max(1,m-n+1)
    i = m
    k = n
    while (i>=l):
        # trouver s'il existe j0 =< k (ligne 4)
        if not all([H[i-1,v]==0 for v in range(k)]):
            # tant que existe j < k avec H[i,j] != 0
            while( not all([H[i-1,v]==0 for v in range(k-1)]) ):
                # trouve ce j0 tant que j<k
                for j in range(1,k):
                    if H[i-1,j-1] != 0:
                        j0 = j
                        break
                # ici on a j0 (ligne 6)
                U.swap_columns(j0-1,k-1)
                H.swap_columns(j0-1,k-1)
                if H[i-1,k-1]<0:
                    H.set_col_to_multiple_of_col(k-1,k-1,-1)
                    U.set_col_to_multiple_of_col(k-1,k-1,-1)
                # pour tous les j < k (ligne 10)
                for jp in range(1,k):
                    c = round(H[i-1,jp-1]/H[i-1,k-1])
                    H.add_multiple_of_column(jp-1,k-1,-c)
                    U.add_multiple_of_column(jp-1,k-1,-c)
            #pour tous les j > k (ligne 12)
            for j in range(k+1,n):
                c = floor(H[i-1,j-1]/H[i-1,k-1])
                H.add_multiple_of_column(j-1,k-1,-c)
                U.add_multiple_of_column(j-1,k-1,-c)
            i = i - 1
            k = k - 1
        else:
            i = i - 1
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
print "La decomposition marche bien, mais elle ne produit le meme résultat que Sage"
print "H == U*A : ", H == A*U

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
    L = identity_matrix(A.base_ring(),m)
    U = identity_matrix(A.base_ring(),n)

    # ECRIVEZ VOTRE CODE ICI
    min_mn = min(m,n)
    # ligne 1
    for k in range(min_mn):
        # ligne 2 -> 7
        if H[k][k] == 0:
            # calcule le produit cartesian
            c = cartesian_product([[k,m-1],[k,n-1]]).list()
            # cherche les indexes tel que Hi0,j0 != 0
            i0 = -1
            j0 = -1
            for v in range (len(c)):
                if (H[c[v][0]][c[v][1]] != 0):
                    i0 = c[v][0]
                    j0 = c[v][1]
                    break
            if i0 != -1:
                H.swap_rows(i0,k)
                L.swap_rows(i0,k)

                H.swap_columns(j0,k)
                U.swap_columns(j0,k)
            else:
                return H
        # ligne 8
        cond = True
        for i in range(k+1,m):
            if H[i][k] != 0:
                cond = False
        for j in range(k+1,n):
            if H[k][j] != 0:
                cond = False
        while(cond == False):
            # check ligne 29
            for i in range (k+1,m):
                xkk = H[k][k]
                xik = H[i][k]
                d,s,t = xgcd(xkk,xik)
                # si xkk divise xik
                if ( xik%xkk == 0):
                    d = xkk
                    s = 1
                    t = 0
                u = -xik/d
                v = xkk/d
                # ligne 17

                Lk = H[k]
                Li = H[i]

                Llk = L[k]
                Lli = L[i]

                # somme terme a terme
                # Lk = s*Lk + t*Li
                H[k] = s*Lk + t*Li
                L[k] = s*Llk + t*Lli
                # Li = u*Lk + v*Li
                H[i] = u*Lk + v*Li
                L[i] = u*Llk + v*Lli
            #ligne 18
            for j in range(k+1,n):
                xkk = H[k][k]
                xkj = H[k][j]
                d,s,t = xgcd(xkk,xkj)
                # si xkk divise xkj
                if ( xkj%xkk == 0):
                    d = xkk
                    s = 1
                    t = 0
                u = -xkj/d
                v = xkk/d
                # les operations suivants sonf faits sur les lignes de la transpose de H
                H = transpose(H)
                U = transpose(U)
                Ck = H[k]
                Cj = H[j]
                Cuk = U[k]
                Cuj = U[j]
                # Ck = s*Ck + t*Cj
                H[k] = s*Ck + t*Cj
                U[k] = s*Cuk + t*Cuj
                # Ci = u*Ck + v*Cj
                H[j] = u*Ck + v*Cj
                U[j] = u*Cuk + v*Cuj
                H = transpose(H)
                U = transpose(U)
            # ligne 27
            # calcule le produit cartesian
            c = cartesian_product([[k,m-1],[k,n-1]]).list()
            # cherche les indexes tel que Hkk |~ Hi0,j0
            i0 = -1
            j0 = -1
            for v in range (len(c)):
                xkk = H[k][k]
                xi0j0 = H[c[v][0]][c[v][1]]
                d,s,t = xgcd(xi0j0,xkk)
                # si reste de xi0j0 par xkk != 0
                if (xi0j0%xkk != 0):
                    j0 = c[v][1]
                    H.add_multiple_of_column(k,j0,1)
                    U.add_multiple_of_column(k,j0,1)
                    break
            #ligne 29
            cond = True
            for i in range(k+1,m):
                if H[i][k] != 0:
                    cond = False
            for j in range(k+1,n):
                if H[k][j] != 0:
                    cond = False
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
print HQ
print "La forme normale de Smith sur F2 est "
print HF2
print "La forme normale de Smith sur F3 est "
print HF3
print "La forme normale de Smith sur F5 est "
print HF5

print "\n$ Question 6"
print "Votre fonction coincide avec celle de Sage ?"
test = True

print '\nX1.smith_form()[0] = \n', X1.smith_form()[0], '\n', 'H1 =\n', H1,'\n'

print '\nX2.smith_form()[0] = \n', X2.smith_form()[0], '\n', 'H2 =\n', H2,'\n'


if AQ.smith_form()[0] != HQ:
    test = False
if AF2.smith_form()[0] != HF2:
    test = False
if AF3.smith_form()[0] != HF3:
    test = False
if AF5.smith_form()[0] != HF5:
    test = False
print test
print 'Again, the result is the same but in some cases it changes the rows or columns.'


reset()
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
H, U = SageHNF(X)
# Solution de Hv = 0
v = matrix(ZZ,[[1],  [0],  [0],  [0]])
# On obtient la representation dans la base originale avec Uv
l = U*v
L =[l] # liste des vecteurs d'une base

# # Affichage des resultats

print "Le systeme a pour racine le module engendre par"
print L

reset()
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
S = A.smith_form()
D = S[0]
L = S[1]
C = S[2]

SageHNF

# dimension of kernel is number of zero columns in H
H, U = SageHNF(A)
H = transpose(H)
k = 0
for i in range(H.nrows()):
    if H[i] == 0:
        k = k + 1

if H.nrows() - k == 4:
    test = True
else:
    test = False

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

D1,L1,C1 = X1.smith_form()

b1_e = L1*b1

z1_e = matrix(ZZ,[
           [ b1_e[0]/D1[0][0]],
           [ b1_e[1]/D1[1][1]],
           [ 0]])

d3_e = matrix(ZZ,[
           [ 0],
           [ 0],
           [ 1]])

d3 = C1*d3_e

z1_ee = C1*z1_e

PolZ.<x_1> = PolynomialRing(ZZ)
z1 = matrix(PolZ,[
           [ z1_ee[0][0] + x_1*d3[0][0]],
           [ z1_ee[1][0] + x_1*d3[1][0]],
           [ z1_ee[2][0] + x_1*d3[2][0]]])

H1 = D1

# Code pour l'EXERCICE

D2,L2,C2 = X2.smith_form()

b2_e = L2*b2

z2_e = D2^(-1)*b2

z2 = C2*z2_e

H2 = D2

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

