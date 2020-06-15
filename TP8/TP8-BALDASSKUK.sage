print("""\
# *************************************************************************** #
# *************************************************************************** #
# TP8 : PRIMALITE DES ENTIERS                                                 #
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
import time

print("""\
# ****************************************************************************
# TEST DE RABIN-MILLER 
# ****************************************************************************
""")

# Code pour l'EXERCICE

# return m and v s.t. n-1 = 2^v*m where m is odd
def get_m_v(n):
    m = 1
    v = 0
    d = n - 1
    while(d%2 == 0):
        d = d/2
        v += 1
    m = d
    return m, v        

def testRM(n):
    a = ZZ.random_element(2,n-1)
    m,v = get_m_v(n)
    g = gcd(a,n)
    if (g > 1):
        return False, g
    b = a^m%n
    if b == 1:
        return True, n
    for i in range(1,v+1):
        if b^2%n == 1:
            g = gcd((b+1), n)
            if g == 1 or g == n:
                return True, n
            else:
                return False, g
        b = b^2%n
    return False, n

def witnessRM(n):
    witnessess = []
    pseudo_primes = []
    m,v = get_m_v(n)
    # test all values in Z/nZ
    for a in range(2,n):
        wit = False
        if mod((a^m),n) != 1:
            wit = True
            for d in range(0,v):
                r = mod(a^((2^d)*m),n)
                if r == n-1:
                    wit = False
                    break
        if wit == True:
            witnessess.append(a)
        if wit == False:
            pseudo_primes.append(a)
    return witnessess, pseudo_primes


for _ in range(10):
    n=ZZ.random_element(3,600)
    print n, "| Resultat Rabin-Miller : ", testRM(n), " | is it truly prime ?", n.is_prime()

print("""\
# ****************************************************************************
#  PERFORMANCES DE RABIN-MILLER
# ****************************************************************************
""")

# Donnees de l'enonce de l'exercice

nmin=10
nmax=500
nbtests = 1

# Code pour l'EXERCICE

n=561

witnessess, pseudo = witnessRM(n)

rep2 = "A witness of 561 is ", witnessess[ZZ.random_element(len(witnessess))], " while it's pseudo-primality bases are ", pseudo
rep3 = "Altough we can clearly see that the excecution time of the algorithm tends to grow, in many cases, the algorithm finds an answer extremely quickly, even for cases where the number is big."


# # Affichage des resultats

bar_chart( [sum( [testRM(n)[1] for i in range(nbtests)]) for n in range(nmin,nmax)]).show()
print rep2
print rep3
list_plot( [timeit( 'testRM(n)', number=20, repeat=3, seconds=true) for n in range(1001,1001+100000,100) ]).show()

reset()
import matplotlib.pyplot as plt

print("""\
# ****************************************************************************
# TEST DE SOLOVAY-STRASSEN 
# ****************************************************************************
""")

# for calculating the Jacobi Symbole
def epsilon(n):
    return (-1)^((n-1)/2)

def omega(n):
    return (-1)^((n^2 - 1)/8)

def upsilon(m,n):
    return (-1)^((m-1)*(n-1)*1/4)

# n >= 3
def Jacobi(m, n):
    if m<0:
        return epsilon(n)*Jacobi(-m,n)
    elif m > 0 and int(m)%2 == 0:
        return omega(n)*Jacobi(m/2,n)
    elif m > 1 and int(m)%2 != 0:
        return upsilon(m,n)*Jacobi(mod(n,m),m)
    elif m == 1:
            return 1
    elif m == 0:
            return 0

def testSS(n):
    if ( n%2 == 0):
        return False
    Zn = Zmod(n)
    a = Zn.random_element()
    while(a <= 1):
        a = Zn.random_element()
    #print'a = '  ,a
    g = gcd(int(a),n)
    if (g != 1):
        return False
    # using sagemath function
    if jacobi_symbol(a,n) == a^((n-1)/2):
        return True
    else:
        return False


rep3 = "A compléter - réponse à la Q3"
rep4 = "A compléter - réponse à la Q4"

# # Affichage des resultats
n = 13
print "Test de la primalite de n= ", n , " avec implementation de Solovay-Strassen ", testSS(n), " | Is n prime: ", n.is_prime()

# # Affichage des resultats
print 'Testing with random integers on range [1000, 1500]'
for _ in range(15):
    n =  ZZ.random_element(1500)
    print 'n is ', n , '| is n prime : ', n.is_prime(), ' | Solovay-Strassen test : ',  testSS(n)

    
print '\nTesting accuracy from range 10, 500 \n'
nmin = 10
nmax = 500
n_tests = 100
N = [i for i in range(nmin, nmax+1)]
T = [0]*len(N)  

k = 0
for i in N:
    #print 'number = ', i+10, '\n'
    for j in range(n_tests):
        #print 'test = ', j
        if ZZ(i).is_prime() == testSS(i):
            T[k] = T[k] + 1
    k += 1
nmin = 10
nmax = 500
n_tests = 100
N = [i for i in range(nmin, nmax+1)]
T = [0]*len(N)  
FP = [0]*len(N) 
k = 0
for i in N:
    #print 'number = ', i+10, '\n'
    for j in range(n_tests):
        #print 'test = ', j
        result = testSS(i)
        if ZZ(i).is_prime() == result:
            T[k] = T[k] + 1
        elif result == True :
            FP[k] += 1
    k += 1
            
fig, ax = plt.subplots()
plt.bar(N, T)
plt.title('Histogram of success')
plt.show()

fig, ax = plt.subplots()
plt.bar(N, FP)
plt.title('Histogram of False Positives')
plt.show()

print 'Some values of false positives are: '
k = 0
while (k < 10):
    n =  ZZ.random_element(len(N))
    if (FP[n] > 0):
        k += 1
        print 'Test failed at ', N[n] , ' prime divisors : ', ZZ(N[n]).prime_divisors()

print '\nSince the test seems to generate false positives with probability < 0.5 for some numbers, we could boost the performance by reaplying the test!\n'


N = [i for i in range(nmin, nmax+1)]
T = [0]*len(N)  
FP = [0]*len(N) 
k = 0
for i in N:
    #print 'number = ', i+10, '\n'
    for j in range(n_tests):
        #print 'test = ', j
        result = (testSS(i) and testSS(i) and testSS(i) and testSS(i))
        if ZZ(i).is_prime() == result:
            T[k] = T[k] + 1
        elif result == True :
            FP[k] += 1
    k += 1
            
if sum(FP) != 0:
    
    fig, ax = plt.subplots()
    plt.bar(N, T)
    plt.title('Histogram of success with test applied 4 times')
    plt.show()
    
    fig, ax = plt.subplots()
    plt.bar(N, FP)
    plt.title('Histogram of False Positives with test applied 4 times')
    plt.show()

    print 'Some values of false positives are: '
    k = 0
    while (k < 10):
        n =  ZZ.random_element(len(N))
        if (FP[n] > 0):
            k += 1
            print 'Test failed at ', N[n] , ' prime divisors : ', ZZ(N[n]).prime_divisors()

else:
    print 'No false positives!'

reset()
print("""\
# ****************************************************************************
# COMPARAISON ENTRE LES TESTS DE R-M ET S-S 
# ****************************************************************************
""")


# Donnees de l'enonce de l'exercice

nmax=150

# Code pour l'EXERCICE

def get_m_v(n):
    m = 1
    v = 0
    d = n - 1
    while(d%2 == 0):
        d = d/2
        v += 1
    m = d
    return m, v  

def witnessSS(n,a):
    js = jacobi_symbol(a,n)
    ap = a^((n-1)/2)
    if mod(js - ap,n) != 0:
        return True
    else :
        return False

def witnessRM(n,a):
    
    m,v = get_m_v(n)
    
    b = mod(a^m,n)
    if b == 1:
        return False
    for i in range(v-1):
        b = mod(b^2,n) 
        if b == -1: 
            return False

    return True

Temoins = []
for n in range(10,151):
    if (n%2 != 0):
        for a in range(n):
            if gcd(a,n)==1 and witnessRM(n,a) and not(witnessSS(n,a)):
                Temoins.append((n,a))
                break

# # Affichage des resultats

print "Liste d'entiers composés et de temoins exclusifs de Rabin-Miller"
print Temoins


reset()
print("""\
# ****************************************************************************
# TEST DE LUCAS
# ****************************************************************************
""")

# returns the result of the sequence uk ad vk using the closed formula
def get_u(k, p, q):
    x = PolynomialRing(RationalField(), 'x').gen()
    f = x^2 - p*x + q
    rho, sigma = f.roots(multiplicities = False)
    return (rho^k - sigma^k)/(rho - sigma)

def get_v(k, p, q):
    x = PolynomialRing(RationalField(), 'x').gen()
    f = x^2 - p*x + q
    rho, sigma = f.roots(multiplicities = False)
    return (rho^k + sigma^k)

# return m and v s.t. n = 2^t*m where m is odd
def get_m_t(c):
    m = 1
    t = 0
    d = c
    while(d%2 == 0):
        d = d/2
        t += 1
    m = d
    return m, t   

# Code pour l'EXERCICE
# returns 0 if composed, 1 if 
def testL(n, p, q):
    delta = p^2 - 4*q
    g = gcd(n,(2*q*delta))
    # 0 : g is a factor of n
    if 1 < g and g < n:
        return False
    # 2 : bad choice of p and q
    elif g == n:
        return False
    
    jac = jacobi_symbol(delta,n)
    m,t = get_m_t(n - jac)
    um = get_u(m, p, q)
    g = gcd(n,um)
    # 0 : g is a factor of n
    if 1 < g and g < n:
        return False
    # 1 : pseudoprime
    elif g == n:
        return True
    for s in range(0,t):
        v2sm = get_v((2^s)*m, p, q)
        g = gcd(n,v2sm)
        if 1 < g and g < n:
            return False
        # 1 : pseudoprime
        elif g == n:
            return True
    # n is composed
    return False
    
# # Affichage des resultats

for _ in range(15):
    n =  ZZ.random_element(50,10000)
    print 'n is ', n , '| is n prime : ', n.is_prime(), ' | Lucas test : ',testL(n, 9, 20)



reset()
print("""\
# ****************************************************************************
# TEST DE BAILLIE, POMERANCE, SELFRIDGE ET WAGSTAFF
# ****************************************************************************
""")

# return m and v s.t. n-1 = 2^v*m where m is odd
def get_m_v(n):
    m = 1
    v = 0
    d = n - 1
    while(d%2 == 0):
        d = d/2
        v += 1
    m = d
    return m, v    

def testRM(n):
    a = 2
    m,v = get_m_v(n)
    g = gcd(a,n)
    if (g > 1):
        return False, g
    b = a^m%n
    if b == 1:
        return True, n
    for i in range(1,v+1):
        if b^2%n == 1:
            g = gcd((b+1), n)
            if g == 1 or g == n:
                return True, n
            else:
                return False, g
        b = b^2%n
    return False, n


# Using the closed formula required complex roots, wich imply on floating point errors.
# for the function to work, the whole sequence of u and v was generated. Not efficient, but easy to implement.
def get_u(m,p,q):
    U = [0]*(int(m)+1)
    U[1] = 1
    k = 2
    while(k <= m):
        U[k] = p*U[k-1] - q*U[k-2]
        k += 1
    return U[int(m)]

def get_v(m,p,q):
    V = [0]*(int(m)+1)
    V[0] = 2
    V[1] = p
    k = 2
    while(k <= m):
        V[k] = p*V[k-1] - q*V[k-2]
        k += 1
    return V[int(m)]


# return m and v s.t. n = 2^t*m where m is odd
def get_m_t(c):
    m = 1
    t = 0
    d = c
    while(d%2 == 0):
        d = d/2
        t += 1
    m = d
    return m, t   

# Code pour l'EXERCICE
# returns 0 if composed, 1 if 
def testL(n, p, q):
    delta = p^2 - 4*q
    g = gcd(n,(2*q*delta))
    # 0 : g is a factor of n
    if 1 < g and g < n:
        return False
    # 2 : bad choice of p and q
    elif g == n:
        return False
    
    jac = jacobi_symbol(delta,n)
    m,t = get_m_t(n - jac)
    um = get_u(m, p, q)
    g = gcd(n,um)
    # 0 : g is a factor of n
    if 1 < g and g < n:
        return False
    # 1 : pseudoprime
    elif g == n:
        return True
    for s in range(0,t):
        v2sm = get_v((2^s)*m, p, q)
        g = gcd(n,v2sm)
        if 1 < g and g < n:
            return False
        # 1 : pseudoprime
        elif g == n:
            return True
    # n is composed
    return False
    

# Code pour l'EXERCICE

def testBPSW(n):
    # avoids a lot of calculations
    if n == 0 or n == 1:
        return False
    if n == 2:
        return True
    elif n%2 == 0:
        return False
    if not testRM(n)[0]:
        return False
    k = 1
    delta = 5
    while jacobi_symbol(delta,n) != -1:
        k+=1
        delta = (-1)^k*(2*k+5)
    p = 1
    q = (1 - delta)/4
    #print 'n = ', n,' p = ', p , ' q = ',q  
    if not testL(n,p,q):
        return False
    return True

nmax=5000

# # Affichage des resultats
print 'Testing with random integers on range [1000, 1500]'
for _ in range(15):
    n =  ZZ.random_element(1500)
    print 'n is ', n , '| is n prime : ', n.is_prime(), ' | BPSW test : ',testBPSW(n)

test = True
for i in range(nmax + 1):
    if ZZ(i).is_prime() != testBPSW(i):
        print '\nn is ', i , '| is n prime : ', ZZ(i).is_prime(), ' | BPSW test : ',testBPSW(i)
        test = False
print '\nDoes the test work for the first ', nmax  ,' integers?\nAnswer: ', test