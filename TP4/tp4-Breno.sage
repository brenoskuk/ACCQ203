print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP4 : FACTORISATION DE POLYNOMES UNIVARIEES SUR CORPS FINIS                 #
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
# FACTORISATION DES PUISSANCES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

F1849.<omega> = FiniteField(43^2,modulus=x^2+1)
Pol1849.<x> = PolynomialRing(F1849)
f=x^172+(3-2*omega)*x^129-5*omega*x^86+(2 + 4*omega)*x^43-1-omega 

F9.<alpha> = FiniteField(9)
Pol9.<y> = PolynomialRing(F9)
g = y^30-y^15+alpha*y^3+1

# Code pour l'EXERCICE

def racine_p_polynome(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    rac = lambda x : x^(q/p)
    u = Pol(map(rac,list(f)[::p]))
    assert(u^p==f)
    return u




# # Affichage des resultats

print "\n$ Question 3"
print "La racine de",g,"est",racine_p_polynome(g)

test = True
# genere 100 polynomes et les teste
for i in range(100):
    p = random_prime(11)
    s = randint(2,5)
    q = p**s
    Fq.<alpha> = FiniteField(q)
    Polp.<y> = PolynomialRing(Fq)
    d_null = False
    # creation du polynome de degree q et puissances multiples de p
    v = [vv for vv in range(q)]
    v = v[0:len(v)-1:p]
    f = 0*y
    for j in v:
        if v != q:
            f = f + Fq.random_element()*y**j
        else: 
            f = f + y**j
        #print 'degree = ',gg.degree(),'\n'
    u = racine_p_polynome(f)
    if (u^p != f):
        print "i = ", i, " ZERO \n"
        test = False
    

print "\n$ Question 4"
print "Test sur 100 exemples : ",test

reset()
print("""\
# ****************************************************************************
# FACTORISATION SANS FACTEURS CARRES
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

F7 = FiniteField(7)
Pol7.<x> = PolynomialRing(F7)
f = x^10 +6*x^9 +3*x^7 +3*x^3 +4*x^2 +2
f2 = 
# Code pour l'EXERCICE



def myFsFC(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    print '\n f = ',f, '\n'
    if (f.degree() <= 0):
        return []
    elif (f.derivative() != 0 ):
        i = 1
        L = []
        t,_,_ = xgcd(f,f.derivative())
        u = f/t
        u = Pol(u)
        while ( u != 1):
            y,_,_ = xgcd(t,u)
            if ((i%p != 0) and (u/y != 1)):
                pp = u/y
                pp = Pol(pp)
                L.append((u/y,i))
            i = i + 1
            u = y
            t = t/y
            t = Pol(t)
        if t != 1:
            r = q/p
            L = L + myFsFC(t^(r))
            
    else:
        L = myFsFC(f^(q/p))
    
    retour = L
    assert(prod([f^e for (f,e) in retour ]) == f)
    return retour

test = false    
    
# # Affichage des resultats

print "\n$ Question 2"
print "La factorisation de",f,"est",myFsFC(f)
print "\n$ Question 4"
print "Test sur 100 exemples : ",test

t = myFsFC(f)


reset()
print("""\
# ****************************************************************************
# FACTORISATION ETAGEE EN DEGRES DISTINCTS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

F5 = FiniteField(5)
Pol.<x>=PolynomialRing(F5)
f = x^10-2*x^9+x^8+x^7-x^6-2*x^5+2*x^4+2*x^3-x

# Code pour l'EXERCICE


def myFEDD(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [f]
    assert(prod(retour) == f)
    return retour

# # Affichage des resultats

print "\n$ Question 1"
print "La factorisation de",f,"est",myFEDD(f)

reset()
print("""\
# ****************************************************************************
# CANTOR-ZASSENHAUSS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

q=3
d=4
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq) 
[f for f in Polq.polynomials(of_degree=d) if f.is_irreducible()
 and f.leading_coefficient()==1]

# Code pour l'EXERCICE

def myCZ(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [f]
    assert(prod(retour) == f)
    return retour
    
def myCZ2(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [f]
    assert(prod(retour) == f)
    return retour    

test  = false

# # Affichage des resultats




reset()
print("""\
# ****************************************************************************
# FACTORISATION COMPLETE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

q=3
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq) 
L1 = [f for f in Polq.polynomials(of_degree=1) if f.is_irreducible()
and f.leading_coefficient()==1]
L2 = [f for f in Polq.polynomials(of_degree=2) if f.is_irreducible()
and f.leading_coefficient()==1]
L3 = [f for f in Polq.polynomials(of_degree=3) if f.is_irreducible()
and f.leading_coefficient()==1]
    
f = L1[0]*L1[1]^3*L1[2]^4
f *= L2[0]*L2[1]^4*L2[2]^4
f *= L3[0]*L3[1]*L3[2]^2*L3[3]^2*L3[4]^3*L3[5]^3*L3[6]^4*L3[7]^4
factor(f)


# Code pour l'EXERCICE

def myFactorisation(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = [(f,1)]
    assert(prod([f^e for (f,e) in retour ]) == f)
    return retour



# # Affichage des resultats

print "\n$ Question 1"
print "La factorisation de",f,"est",myFactorisation(f)



reset()
print("""\
# ****************************************************************************
# RACINES D'UN POLYNOME
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

q=3
Fq=FiniteField(q)
Polq.<x> = PolynomialRing(Fq) 
L1 = [f for f in Polq.polynomials(of_degree=1) if f.is_irreducible()
and f.leading_coefficient()==1]
L2 = [f for f in Polq.polynomials(of_degree=2) if f.is_irreducible()
and f.leading_coefficient()==1]
L3 = [f for f in Polq.polynomials(of_degree=3) if f.is_irreducible()
and f.leading_coefficient()==1]
    
f = L1[0]*L1[1]^3*L1[2]^4
f *= L2[0]*L2[1]^4*L2[2]^4
f *= L3[0]*L3[1]*L3[2]^2*L3[3]^2*L3[4]^3*L3[5]^3*L3[6]^4*L3[7]^4


# Code pour l'EXERCICE


def myRacine(f):
    Pol=f.parent()
    x=Pol.gen()
    p=Pol.base_ring().characteristic()
    q=Pol.base_ring().cardinality()
    retour = []
    assert(f(z)==0 for z in retour)
    return retour

# # Affichage des resultats

print "\n$ Question 1"
print "Les racines de ",f,"sont",myRacine(f)




reset()
print("""\
# ****************************************************************************
# ETUDE DE CANTOR-ZASSENHAUSS
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

# Code pour l'EXERCICE


# # Affichage des resultats

