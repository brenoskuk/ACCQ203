print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP13 : COURBES ELLIPTIQUES                                                  #
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
# ADDITION DANS UNE COURBE ELLIPTIQUE
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 61
Fp = FiniteField(p)
E = EllipticCurve(Fp,[1,2])
P = E.random_point()
Q = E.random_point()

# Code pour l'EXERCICE


def addition(P,Q):
    
    E = P.curve()
    Fp = E.base_ring()
    
    if Q.is_zero():
        return P    
    else:
        xq,yq = Q.xy()
        
    if P.is_zero():
        return Q
    else:
        xp,yp = P.xy()
    
    # different formulas
    if P != Q:     
        # point at infinite
        if xp == xq:
            return E.point([0,1,0])

        lamb = (yq-yp)/(xq-xp)
        ni = (yp*xq - yq*xp)/(xq - xp)
        a = lamb^2 - xp - xq
        b = -lamb^3 + lamb*(xp + xq) - ni
        return E.point([a,b])

    else:

        if yq != 0:
            lamb = (3*xp^2 + E.a4())/(2*yp)
            ni = (-3*x^3 + E.a4()*x + 2*E.a6())/(2*yp)
            
            a = lamb^2 - xp - xq
            b = -lamb^3 + lamb*(xp + xq) - ni

            return E.point([a,b])
        
        else:
            return E.point([0,1,0])

    
    
# # Affichage des resultats
for _ in range(10):
    print 'P = ', P, ' Q = ', Q 
    print 'My addition: P + Q = ', addition(P,Q), '\nSagemath :   P + Q = ', P + Q, '\n'


reset()
print("""\
# ****************************************************************************
# COURBE DE L'ANSSI
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE

ANSSI = "Agence nationale de la sécurité des systèmes d'information"
p = 109454571331697278617670725030735128145969349647868738157201323556196022393859
a = 82638672503301278923015998535776227331280144783487139112686874194432446389503
b = 43992510890276411535679659957604584722077886330284298232193264058442323471611
E = EllipticCurve(FiniteField(p),[a,b])




# # Affichage des resultats

print "ANSSI signifie :",ANSSI
print "La courbe recommandée est"
print E 
print p.is_prime()
print "\nG = ", E.cardinality() ," and log2(G) = " , float(log(sqrt(E.cardinality()),2)), " " 

reset()
print("""\
# ****************************************************************************
# COMPTAGE DE POINTS
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

p = 2003
Fp = FiniteField(p)
a = 1929
b = 1178

while true:
    d=Fp.random_element()
    if not d.is_square():
        break


# Code pour l'EXERCICE


def comptage(a,b,Fp):  
    # at least 0e
    count = 1
    # test all elements of Fp
    for x in Fp:
        y2 = x^3 + a*x+ b
        if y2.is_square():
            count = count + 2
        elif y2 == 0:
            count = count + 1
    
    return count

frequence = [4,1,4]

# # Affichage des resultats

E = EllipticCurve(Fp, [a,b])
print "My count = ", comptage(a,b,Fp), " Sagemath = ", E.cardinality()

# for a finite field Fq returns freq of orders
def histo_orders(q):
    Fq = FiniteField(q)
    # get max order and save orders if they exist:
    max_order = 0
    orders = []
    for a in Fq:
        for b in Fq:
            try:
                order = EllipticCurve(Fq,[a,b]).cardinality()
            except:
                order = -1
            if order > max_order:
                max_order = order
            if order > 0:
                orders.append(order)

    freq = [0]*max_order

    for order in orders:
        freq[order-1]+=1

    return freq
    
q = random_prime(200)
print "For q = ", q
print "Hasse's threshold: ", q + 1 + int(sqrt(p))
chart = bar_chart(histo_orders(q))
chart.plot()
chart.show()

reset()
print("""\
# ****************************************************************************
# FACTORISATION ECM
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

class FoundFactor(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

n = 2020

# Code pour l'EXERCICE

def division(x,y):
    try:
        r = x/y
        return r
    except:
        if y == 0:
            raise ZeroDivisionError
        else:
            raise FoundFactor(y)
    
    
    
def addition(P,Q):
    
    E = P.curve()
    Fp = E.base_ring()
    
    if Q.is_zero():
        return P    
    else:
        xq,yq = Q.xy()
        
    if P.is_zero():
        return Q
    else:
        xp,yp = P.xy()
    
    # different formulas
    if P != Q:     
        # point at infinite
        if xp == xq:
            return E.point([0,1,0])

        lamb = division((yq-yp),(xq-xp))
        ni = (yp*xq - yq*xp)/(xq - xp)
        a = lamb^2 - xp - xq
        b = -lamb^3 + lamb*(xp + xq) - ni
        return E.point([a,b])

    else:

        if yq != 0:
            lamb = division((3*xp^2 + E.a4()),(2*yp))
            ni = (-3*x^3 + E.a4()*x + 2*E.a6())/(2*yp)
            
            a = lamb^2 - xp - xq
            b = -lamb^3 + lamb*(xp + xq) - ni

            return E.point([a,b])
        
        else:
            return E.point([0,1,0])
        

def multiplication(lamb,P):
    
    if lamb == 0:
        return P.parent()(0)
    if lamb == 1:
        return P
    if lamb % 2 == 0:
        return multiplication((lamb/2),addition(P,P))
    if lamb % 2 == 1:
        return addition(P , multiplication((lamb - 1)/2,addition(P,P)))
        


def get_e(p,B):
    e = 1
    while (p^e < B):
        e += 1
    return e

def ECM(n,B):
    
    Zn = Integers(n)
    
    a = Zn.random_element()
    x0 = Zn.random_element()
    y0 = Zn.random_element()
    
    b = y0^2 - x0^3 - a*x0
    
    g = gcd(4*a^3+27*b^3 , n)
    
    if (1<g) and (g < n):
        return g
    elif g == n :
        return False
    
    E = EllipticCurve(Zn,[a,b])
    A = E.point([x0,y0])
    
    P = primes(B)
    for p in P:
        e = get_e(p,B)
        
        try:
            Ap = multiplication(p^e, A)
            A = Ap
        except FoundFactor as FF:
            return gcd(FF.value,n)
            
            
    return False
    


# # Affichage des resultats

n = 2020
B = 500
print 'n = ', n
res = ECM(n,B)
print 'ECM(n,500) = ', res
print 'division result = ', n/int(res), '\n'
i = 0
e = 0
while i < 10:
    n = randint(100,4000)
    try:
        res = ECM(n,B)
        i = i + 1
        print 'n = ', n
        if res == False:
            'Echoue'
        else:
            print 'ECM(n,500) = ', res
            print 'division result = ', n/int(res), '\n'
    except:
        e += 1


reset()
print("""\
# ****************************************************************************
# EXPONENTIATION RAPIDE ET ATTAQUE PAR CANAUX CACHES
# ****************************************************************************
""")

# NE PAS TRAITER



reset()
print("""\
# ****************************************************************************
# COURBE D'EDWARDS
# ****************************************************************************
""")



# Code pour l'EXERCICE

reponse = "With an Edwards curve, we calculate a duplication as we would calculate an addition. Thus we can't infer the operations using their resulting noise."

# # Affichage des resultats

print reponse