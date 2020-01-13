print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP6 : BASES DE GROEBNER ET SYSTEMES POLYNOMIAUX MULTIVARIES                 #
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
#  FONCTIONS DE SAGEMATH
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x,y,z> = PolynomialRing(QQ,3, order='lex')
f = 2*x^2*y+7*z^3

# Code pour l'EXERCICE

print x<y^2
print f.lt()
print f.lc()
print f.lm()

reponse  ="\n The '<' comparison is a monomial order (defined as lexicographical)\n f.lt() is the dominant term\n f.lc() is the dominant coefficient\n f.lm() is the dominant monimial\n The documantation cites the following orders(besides lexicographical): \n  -Degree reverse lexicographic (degrevlex)\n  -Degree lexicographic (deglex)\n  -Inverse lexicographic (invlex)\n  -Negative lexicographic (neglex)\n  -Negative degree reverse lexicographic (negdegrevlex)\n  -Negative degree lexicographic (negdeglex)\n  -Weighted degree reverse lexicographic (wdegrevlex), positive integral weights\n  -Weighted degree lexicographic (wdeglex), positive integral weights\n  -Negative weighted degree reverse lexicographic (negwdegrevlex), positive integral weights\n  -Degree negative lexicographic (degneglex)\n  -Negative weighted degree lexicographic (negwdeglex), positive integral weights"

# # Affichage des resultats

print "\n$1/ ", reponse

reset()
print("""\
# ****************************************************************************
# DIVISION MULTIVARIEE
# ****************************************************************************
""")




# Donnees de l'enonce de l'exercice

MPol.<x,y> = PolynomialRing(QQ,2, order='lex')
f  = -x^7 + x^6*y + 2*x^5 - 2*x^4*y - 5*x^2 + 3*x*y^3 + 5*x*y + 11*y^3 + 10 
f1 = x*y^2+2*y^2
f2 = x^5+5

# Code pour l'EXERCICE

def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    Q = [MPol(0)]*s
    r = f
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r

# # Affichage des resultats

print "\n$ ",  myDivision(f,[f1,f2])

reset()
print("""\
# ****************************************************************************
# BASE DE GROEBNER
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

MPol.<x,y,z> = PolynomialRing(QQ,3, order='lex')
f1 = x^2-y
f2 = x*y-z
f3 = z^4+x*y

# Code pour l'EXERCICE
def Syzygie(g,h):
    MPol = g.parent()
    tdg = g.lt()
    tdh = h.lt()
    S = lcm(tdg,tdh)*g/tdg - lcm(tdg,tdh)*h/tdh
    S = MPol(S)
    return S

def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    # line 1
    r = 0
    p = f
    # line 2
    Q = [MPol(0)]*s
    # line 4
    while (p != 0):
        # line 5
        cond = False 
        for i in range(0,s): 
            if (p.lt()%F[i].lt() == 0):
                Q[i] = Q[i] + p.lt()/F[i].lt()
                p = p - MPol((p.lt()/F[i].lt())*F[i])
                cond = True
                break
        # line 8
        if not(cond):
            r = r + p.lt()
            p = p - p.lt()
        #print 'p = ',p,'\n','r = ', r, '\nQ = ', Q, '\np.lt()', p.lt(),'\n'
        
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r

# Algo 16
def myGroebner(F):
    # line 1
    G = copy(F)
    S_void = False
    # line 4
    while not(S_void):
        # line 3
        S=[]
        S_void = True
        # line 4
        for i in range(len(G)):
            for j in range(i+1,len(G)):
                # line 5
                Syz = Syzygie(G[i],G[j])
                r = myDivision(Syz,G)[1]
                # line 6 and 7
                if r != 0:
                    S = list(set(S) | set([r])) 
                    S_void = False
        G = list(set(G) | set(S))
        G.sort()
    return G
    
def myRedGroebner(G):
    i = 0
    while (i < len(G)):
        G_t = copy(G)
        g = G_t.pop(i)
        r = myDivision(g,G_t)[1]
        if (r == 0):
            G = G_t
        else:
            i+=1
    # normalisation
    for i in range(len(G)):
        G[i] = G[i]/(G[i].coefficients()[0])
    return G

# # Affichage des resultats
G0 = myGroebner([f1,f2,f3])
print "\n$1/ ",G0
Gr = myRedGroebner(G0)
print "\n$2/ ",Gr



reset()
print("""\
# ****************************************************************************
# APPARTENANCE A UN IDEAL
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

MPol.<x,y,z> = PolynomialRing(QQ,3, order='lex')
f1 = x*y-y^2
f2 = x^3-z^2
I = Ideal([f1,f2])
f = -4*x^2*y^2*z^2 + y^6 + 3*z^5

# Code pour l'EXERCICE

def Syzygie(g,h):
    MPol = g.parent()
    tdg = g.lt()
    tdh = h.lt()
    S = lcm(tdg,tdh)*g/tdg - lcm(tdg,tdh)*h/tdh
    S = MPol(S)
    return S

def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    # line 1
    r = 0
    p = f
    # line 2
    Q = [MPol(0)]*s
    # line 4
    while (p != 0):
        # line 5
        cond = False 
        for i in range(0,s): 
            if (p.lt()%F[i].lt() == 0):
                Q[i] = Q[i] + p.lt()/F[i].lt()
                p = p - MPol((p.lt()/F[i].lt())*F[i])
                cond = True
                break
        # line 8
        if not(cond):
            r = r + p.lt()
            p = p - p.lt()
        #print 'p = ',p,'\n','r = ', r, '\nQ = ', Q, '\np.lt()', p.lt(),'\n'
        
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r

# Algo 16
def myGroebner(F):
    # line 1
    G = copy(F)
    S_void = False
    # line 4
    while not(S_void):
        # line 3
        S=[]
        S_void = True
        # line 4
        for i in range(len(G)):
            for j in range(i+1,len(G)):
                # line 5
                Syz = Syzygie(G[i],G[j])
                r = myDivision(Syz,G)[1]
                # line 6 and 7
                if r != 0:
                    S = list(set(S) | set([r])) 
                    S_void = False
        G = list(set(G) | set(S))
        G.sort()
    return G
    
def myRedGroebner(G):
    i = 0
    while (i < len(G)):
        G_t = copy(G)
        g = G_t.pop(i)
        r = myDivision(g,G_t)[1]
        if (r == 0):
            G = G_t
        else:
            i+=1
    # normalisation
    for i in range(len(G)):
        G[i] = G[i]/(G[i].coefficients()[0])
    return G

G = myRedGroebner(myGroebner([f1,f2]))

_,rr = myDivision(f,G)

if rr == 0:
    test1 = True
else:
    test1 = False
    
test1 = f in I
n = randint(1,10)
ff = I.random_element(degree=n)

_,rr = myDivision(ff,G)


if rr == 0:
    test2 = True
else:
    test2 = False


# # Affichage des resultats

print "\n$ Test de Sage ",test1
print "\n$ Test de personnel ",test2


reset()
print("""\
# ****************************************************************************
# RESOLUTION D'UN SYSTEME
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


MPol.<x,y> = PolynomialRing(QQ,2) # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
f = (y^2+6)*(x-1) - y*(x^2 + 1)
g = (x^2+6)*(y-1) - x*(y^2 + 1)
 

# Code pour l'EXERCICE

def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    # line 1
    r = 0
    p = f
    # line 2
    Q = [MPol(0)]*s
    # line 4
    while (p != 0):
        # line 5
        cond = False 
        for i in range(0,s): 
            if (p.lt()%F[i].lt() == 0):
                Q[i] = Q[i] + p.lt()/F[i].lt()
                p = p - MPol((p.lt()/F[i].lt())*F[i])
                cond = True
                break
        # line 8
        if not(cond):
            r = r + p.lt()
            p = p - p.lt()
        #print 'p = ',p,'\n','r = ', r, '\nQ = ', Q, '\np.lt()', p.lt(),'\n'
        
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r

base = Ideal([f,g]).groebner_basis()

# est la est reduite? 
_,rr = myDivision(base[1],[base[0],base[1]])

if rr == 0:
    base = [base[0],base[2]]


py0 = base[0].univariate_polynomial()

racines_y = py0.roots(multiplicities=False) 

# ici on a un sous ensemble de racines de y

px0_1 = base[1].subs(y = racines_y[0]).univariate_polynomial()
px0_2 = base[1].subs(y = racines_y[1]).univariate_polynomial()


# px0_1 == px0_2 est vrai, donc on travaile juste avec px0_1

racines_x = px0_1.roots(multiplicities=False)

racines = cartesian_product([racines_x,racines_y]).list()

Gf = implicit_plot(f,(x,0,6),(y,0,6),color='red') 
Gg = implicit_plot(g,(x,0,6),(y,0,6),color='blue')  
Gp = point2d(racines,color='green')

# # Affichage des resultats

print "\n$1/  Une base de Groebner de [f,g] est", base
print "\n$2/  Les valeurs de y sont", racines_y
print "\n$4/  Les valeurs de (x,y) sont", racines
print "\n$4/"; show(Gf+Gg+Gp)

reset()
print("""\
# ****************************************************************************
# OPTIMISATION
# ****************************************************************************
""")


# Donnees et code de l'exercice


# selon la Proposition 254. on choisit l'ordre lexicographique

MPol.<x,y,lamb> = PolynomialRing(QQ,3,order='lex') # QUEL ORDRE DEVEZ-VOUS CHOISIR ?
f = x^2*y  - 2*x*y + y + 1
g = x^2 + y^2 - 1

# on calcule les gradients

gradf = f.gradient()
gradg = g.gradient()

# on ajoute la condition de Karash Khun Tucker en chaque composant

kx = gradf[0] - lamb*gradg[0]
ky = gradf[1] - lamb*gradg[1]

# on va reçoudre le systeme multivariable:
# kx = 0
# ky = 0
# g = 0
# dans l'ideal avec Grobner

syst = [str(g) + ' = 0',str(kx) + ' = 0',str(ky) + ' = 0']

base = Ideal([g,kx,ky]).groebner_basis()

# le dernier polynome de la base est seulement fonction de lambda

plamb = base[-1].univariate_polynomial()

# on cherche les racines de lambda sur RR
racines_plamb = plamb.roots(ring = RR, multiplicities=False) 

racines = []
for rac in racines_plamb:
    L1 = (base[1].subs(lamb=rac)).univariate_polynomial().roots(ring=RR, multiplicities=false) 
    for y0 in L1:
        for x0 in (base[0].subs(lamb=rac, y =y0)).univariate_polynomial().roots(ring=RR, multiplicities=false):
            racines.append((x0,y0))


# on cherche le maxima dans f dans les racines
sol = racines[0]
fmax = f.subs(x = sol[0], y = sol[1])
for rac in racines:
    fp = f.subs(x = rac[0], y = rac[1])
    if fp>fmax:
        fmax = fp
        sol = rac

# # Affichage des resultats


print "\n$1/  On doit resoudre le systeme\n", 
for i in range(len(syst)):
    print syst[i]
print "\n$2/  dont une base de Groebner est\n",
for i in range(len(base)):
    print '[',base[i],']'
print "\n$4/  Les valeurs de (x,y) sont", sol


xx, yy = var('x y')
from sage.plot.plot3d.plot3d import axes
W = plot3d(xx^2*yy - 2*xx*yy^2 + yy + 1,(x,-1.5,1.5),(y,-1.5,1.5), frame=False, color='purple', opacity=0.8)
S=Cylindrical('radius', ['azimuth', 'height'])
theta,z=var('theta, z')
C = plot3d(1, (theta,0,2*pi), (z, -2, 3), transformation=S)
M = sphere((sol[0],sol[1],fmax), size=.2, color='red') 
show(W + C + M + axes(2, color='black'), figsize=8)

reset()
print("""\
# ****************************************************************************
# MANIPULATIONS ALGEBRIQUES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

# Code pour l'EXERCICE

def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    # line 1
    r = 0
    p = f
    # line 2
    Q = [MPol(0)]*s
    # line 4
    while (p != 0):
        # line 5
        cond = False 
        for i in range(0,s): 
            if (p.lt()%F[i].lt() == 0):
                Q[i] = Q[i] + p.lt()/F[i].lt()
                p = p - MPol((p.lt()/F[i].lt())*F[i])
                cond = True
                break
        # line 8
        if not(cond):
            r = r + p.lt()
            p = p - p.lt()
        #print 'p = ',p,'\n','r = ', r, '\nQ = ', Q, '\np.lt()', p.lt(),'\n'
        
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r

MPol.<c,s,u,v> = PolynomialRing(QQ,4,order='lex') 

# s = sin ; c = cos -> v = 2sc + c^2 - s^2; s^2 + c^2 = 1

f1 = u - s - c
f2 = v - 2*s*c - c^2 + s^2
f3 = s^2 + c^2 - 1

# A partir de l'ideal R = [f1,f2,f3] on trouve une base de Groebner

G = Ideal([f1,f2,f3]).groebner_basis()

# Finalmente on projete s^6 sur cette base

_,r = myDivision(s^6,G)

# on note que pour n'importe quel theta, si p est dans l'ideal R
# p*f1 = p*f2 = p*f3 = 0
# une bonne donc aproximation est donne par le reste
# qui n'aperient pas a R et consist seulement de u e v:

# # Affichage des resultats

print "\n$1/ "

print "s^6 ~= ",r

# test
#tt = 1
#cc = float(cos(tt))
#ss = float(sin(tt))
#ss6 = ss^6
#uu = ss + cc
#vv = 2*ss*cc + cc^2 - ss^2
#q4 = 1/8*(uu^2) - 1/4*(vv) + 1/4
#rr = 1/16*uu^2*vv^2 - 3/8*uu^2*vv + 7/16*uu^2 + 1/8*vv^2 - 1/8*vv - 1/8

reset()
print("""\
# ****************************************************************************
# OVALES DE DESCARTES
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice


# Code pour l'EXERCICE

def myDivision(f,F):
    MPol = f.parent()
    n = MPol.ngens()
    s = len(F)
    # line 1
    r = 0
    p = f
    # line 2
    Q = [MPol(0)]*s
    # line 4
    while (p != 0):
        # line 5
        cond = False 
        for i in range(0,s): 
            if (p.lt()%F[i].lt() == 0):
                Q[i] = Q[i] + p.lt()/F[i].lt()
                p = p - MPol((p.lt()/F[i].lt())*F[i])
                cond = True
                break
        # line 8
        if not(cond):
            r = r + p.lt()
            p = p - p.lt()
        #print 'p = ',p,'\n','r = ', r, '\nQ = ', Q, '\np.lt()', p.lt(),'\n'
        
    assert(f==sum(q*g for q,g in zip(Q,F) )+r)
    return Q,r


MPol.<OM,BM,x,y> = PolynomialRing(QQ,4,order='lex') 

# on ecrit les equations
f1 = OM^2 - x^2 - y^2
f2 = BM^2 - (x-1)^2 + y^2
f3 = OM + 2*BM - 3
 
# On sait que les solutions du sisteme ci dessus sont aussi
# de la base de Groebner
G = Ideal([f1,f2,f3]).groebner_basis()


# Il y a une solution dans G qui contient seulement des coordones x,y
# cela c'est lieu des points M cherché

eq = G[2]

# Affichage des resultats
MPlot.<x,y> = PolynomialRing(QQ,2) 
eq = MPlot(eq)
Gf = implicit_plot(eq,(x,-10,10),(y,-10,10),color='red') 
print "\n$4/"; show(Gf)



print "\n$ L'équation est ",eq, ' = 0'

