print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP3 : RESEAUX EUCLIDIENS                                                    #
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
# BASE D'UN RESEAU
# ****************************************************************************
""")


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

# Donnees de l'enonce de l'exercice

A = matrix(ZZ,[
        [ 2, -3,  4, 0],
        [ 4,  4, 53, -5]])

# Code pour l'EXERCICE

H,U = SageHNF(A)


#Solutions of Ax = b -> Ax


print H, '\n', U, '\n', H - A*U == 0

# prop 57 - nombre de colomnes zero
L = []
L.append(vector(U[:3,0]))
L.append(vector(U[:3,1]))

# # Affichage des resultats

print "\n$ Le réseau a pour base"
print L



reset()
print("""\
# ****************************************************************************
# APPLICATIONS NUMERIQUES
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

n1 = round(arctan(1),50)
n2 = round(arctan(1/5),50)
n3 = round(arctan(1/239),50)

r=-2.5468182768840820791359975088097915

Pol.<x>=PolynomialRing(ZZ)

# Code pour l'EXERCICE
M = 1000000
B = matrix(ZZ,[[1,0,0],
              [0,1,0],
              [0,0,1],
               [round(M*n1),round(M*n2),round(M*n3)]])

B = transpose(B)
alpha = B.LLL()
alpha = transpose(alpha)

B2 = matrix(ZZ,[[1,0,0,0],
              [0,1,0,0],
              [0,0,1,0],
              [0,0,0,1],
               [round(M),round(M*r),round(M*r^2),round(M*r^3)]])

B2 = transpose(B2)
alpha2 = B2.LLL()
alpha2 = transpose(alpha2)
Pol.<x>=PolynomialRing(QQ)
p = alpha2[0,0] + alpha2[1,0]*x + alpha2[2,0]*x^2 + alpha2[3,0]*x^3
rr = p.roots()


# # Affichage des resultats

print "\n$ La relation de Machin est alpha1*n1+alpha2*n2+alpha3*n3=0 avec"
for i in range(3):
   print "alpha",i+1,"=",alpha[i]

print "\n$ Un polynome minimal plausible est"
print p
print "dont les racines sont"
print p.roots(ring=RR,multiplicities=false)




reset()
print("""\
# ****************************************************************************
# ALGORITHME LLL
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

B = matrix(ZZ,[[  9, -25],
               [ -5,  14],
               [  4, -11]])

# Code pour l'EXERCICE

# Sans mise a jour
def myLLL(M):
    # ligne 1
    repeter = True
    B = copy(M)
    n = B.ncols()
    while (repeter):
        # ligne 3
        GH, mu = B.transpose().gram_schmidt(orthonormal=False)
        # On ne besoin de GH dans la convention française
        # GH = GH.transpose()
        mu = mu.transpose()
        # ligne 4
        for i in range(1,n):
            for j in range(i-1,-1,-1):
                # transposition pour faire des operations en ligne
                mji = mu[j,i]
                B = B.transpose()
                mu = mu.transpose()
                B[i] = B[i] - round(mji)*B[j]
                mu[i] = mu[i] - round(mji)*mu[j]
                # transposition pour retourner a B originale
                mu = mu.transpose()
                B = B.transpose()
        # ligne 8
        echange = False
        i0 = 0
        while((i0 < GH.nrows() - 1) and not(echange)):
            if (abs(GH[i0])^2 > 2*abs(GH[i0+1])^2):
                B = B.transpose()
                L = B[i0]
                B[i0] = B[i0+1]
                B[i0+1] = L
                B = B.transpose()
                echange = True
            i0 = i0 + 1
        if (echange == False):
            repeter = False
    return B

# Avec mise a jour
def myLLL_ms(M):
    # ligne 1
    repeter = True
    B = copy(M)
    n = B.ncols()
    GH, mu = B.transpose().gram_schmidt(orthonormal=False)
    # On ne besoin de GH dans la convention française
    # GH = GH.transpose()
    mu = mu.transpose()
    while (repeter):
        for i in range(1,n):
            for j in range(i-1,-1,-1):
                # transposition pour faire des operations en ligne
                mji = mu[j,i]
                B = B.transpose()
                mu = mu.transpose()
                B[i] = B[i] - round(mji)*B[j]
                mu[i] = mu[i] - round(mji)*mu[j]
                # transposition pour retourner a B originale
                mu = mu.transpose()
                B = B.transpose()
        # ligne 7
        echange = False
        i0 = 0
        while((i0 < GH.nrows() - 1) and not(echange)):
            if (abs(GH[i0])^2 > 2*abs(GH[i0+1])^2):
                echange = True
                B = B.transpose()
                L = B[i0]
                B[i0] = B[i0+1]
                B[i0+1] = L
                B = B.transpose()
                s = GH[i0]
                t = GH[i0+1] + mu[i0,i0+1]*GH[i0]
                # ligne 13
                mu = mu.transpose()
                Lmu = B[i0]
                mu[i0] = mu[i0+1]
                mu[i0+1] = Lmu
                mu = mu.transpose()
                # ligne 14
                GH[i0] = t
                muii1 = s.dot_product(GH[i0])/GH[i0].dot_product(GH[i0])
                mu[i0,i0+1] = muii1
                GH[i0] = s - muii1*GH[i0]
                for k in range(i0+2,n):
                    bkbi = GH[k].dot_product(GH[i0])/GH[i0].dot_product(GH[i0])
                    bkbi1 = GH[k].dot_product(GH[i0+1])/GH[i0].dot_product(GH[i0+1])
                    mu[i0,k] = bkbi
                    mu[i0+1,k] = bkbi1
            i0 = i0 + 1
            if (echange == False):
                repeter = False
    return B

# # Affichage des resultats

print "\n$ Une base LLL de B est"
print myLLL(B)


reset()
print("""\
# ****************************************************************************
# RESEAUX CLASSIQUE
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

n = 8
e = vector([1/2]*8)

# Code pour l'EXERCICE

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


# 1.
An = matrix(ZZ,[1,1,1,1,1,1,1,1,1])
H1,U1 = SageHNF(An)

# prop 57 - nombre de colomnes zero
r1 = 0
while (H1.transpose()[r1] == 0):
    r1 = r1 + 1

K1 = U1[:,:r1]
# det of basis is det(Bt*B)
an = det(K1.transpose()*K1) #determinant de An

# 2.
Dn = matrix(ZZ,[1,1,1,1,1,1,1,1,-2])
H2,U2 = SageHNF(Dn)

# prop 57 - nombre de colomnes zero
r2 = 0
while (H2.transpose()[r2] == 0):
    r2 = r2 + 1

K2 = U2[:,:r2]
# det of basis is det(Bt*B)
dn = det(K2.transpose()*K2) #determinant de An

# 3.
E8 = matrix(QQ,[[2,-1,-1,-1,-1,-1,-1,1/2],
                [0,1,0,0,0,0,0,1/2],
                [0,0,1,0,0,0,0,1/2],
                [0,0,0,1,0,0,0,1/2],
                [0,0,0,0,1,0,0,1/2],
                [0,0,0,0,0,1,0,1/2],
                [0,0,0,0,0,0,1,1/2],
                [0,0,0,0,0,0,0,1/2]])

H3,U3 = SageHNF(E8)
# prop 57 - nombre de colomnes zero
r3 = 0
while (H3.transpose()[r3] == 0):
    r3 = r3 + 1

K3 = U3[:,:r3]
# det of basis is det(Bt*B)
e8 = det(K3.transpose()*K3) #determinant de An

reponse6 = "Toutes les elements de la base ont des '1' apart de un qui a seulement de '1/2'. Si il est non nule le nombre resultant es demi-entier. Sinon il est entier "
reponse7_an = "Toutes vecteurs qui ont deux positions +-1 et le reste zero"
reponse7_dn = ""
reponse7_e8 = "Union du cas pour "
reponse8 = "les vecteurs minimaux sont les vecteurs de norme carre 2"

# # Affichage des resultats

print "\n$ Une base de An est"
print K1, "de déterminant",an

print "\n$ Une base de Dn est"
print K2, "de déterminant",dn

print "\n$ Une base de E8 est"
print E8, "de déterminant",e8


reset()
print("""\
# ****************************************************************************
# DENSITES OPTIMALES
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE

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


O4 = matrix(ZZ,[1,1,1,1,-2])
D4,C4 = SageHNF(O4)
r4 = 0
while (D4.transpose()[r4] == 0):
    r4 = r4 + 1

K4 = C4[:,:r4]

det4 = det(K4.transpose()*K4)

O5 = matrix(ZZ,[1,1,1,1,1,-2])
D5,C5 = SageHNF(O5)
r5 = 0
while (D5.transpose()[r5] == 0):
    r5 = r5 + 1

K5 = C5[:,:r5]

det5 = det(K5.transpose()*K5)

# directement de la definition et en prenant les vecteur minimaux
# les vecteurs minimaux sont les vecteurs de norme carre 2

def d(n):
    return pi^(n/2)/(gamma(n/2)*2**n)*sqrt(2)**n


a2 = d(2)/(3**0.5)
a3 = d(3)/(4**0.5)
d4 = d(4)/(det4**0.5)
d5 = d(5)/(det5**0.5)
e8 = d(8)

# # Affichage des resultats

print "\n$ La densité de A2 est",a2
print "\n$ La densité de A3 est",a3
print "\n$ La densité de D4 est",d4
print "\n$ La densité de D5 est",d5
print "\n$ La densité de E8 est",e8




