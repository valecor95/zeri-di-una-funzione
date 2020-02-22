import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


''' ------------------------------------------------- METODO DELLA BISEZIONE '''

def bisezione(f,a,b,nmax,toll):
    xvect = np.array([])
    if f(a)*f(b) >= 0:
        print("ERRORE: f(a) and f(b) devono avere valori differenti")
        return None
    a_n = a
    b_n = b
    it=0
    while(it<nmax and abs(b_n-a_n)>toll):
        c = (a_n + b_n)/2
        xvect = np.append(xvect, [c])
        fc = f(c)
        if f(a_n)*fc < 0:
            a_n = a_n
            b_n = c
        elif f(b_n)*fc < 0:
            a_n = c
            b_n = b_n
        elif fc == 0:
            print("Trovata la soluzione esatta")
            return c
        else:
            print("Errore.")
            return None
        it+=1
    return xvect


''' ------------------------------------------------------- METODO DI NEWTON '''

def newton(f,Df,x0,nmax,toll):
    xvect = np.array([])
    xn = x0
    it = 0

    while(it<nmax):
        fxn = f(xn)
        if abs(fxn) < toll:
            print('Soluzione trovata dopo',it,'iterazioni.')
            return xvect
        Dfxn = Df(xn)
        if Dfxn == 0:
            print('Derivata uguale a 0. Nessuna soluzione trovata.')
            return None
        xn = xn - fxn/Dfxn
        xvect = np.append(xvect, [xn])
        it += 1
    print('Numero di massimo di iterazioni superato. Nessuna soluzione trovata.')
    return None


''' --------------------------------------------------- METODO DELLA SECANTE '''

def secante(f,a,b,nmax,toll):
    xvect = np.array([])
    a_n = a
    b_n = b
    it = 0
    while(it<nmax and abs(b_n-a_n)>toll):
        fb = f(b_n)
        fa = f(a_n)
        c = (fb-fa)/(b_n-a_n);
        if(fb == 0):
            break
        x2 = b_n - (fb/c)
        a_n = b_n
        b_n = x2
        xvect = np.append(xvect, [b_n])
        it += 1
    return xvect



''' ------------------------------------------------------------------- TEST '''

f = lambda x: x**2 - x - 1                                                      # Funzione di input
Df = lambda x: 2*x - 1                                                          # Derivata della funzione di input

print("#################################### BISEZIONE")
my_bis_res = bisezione(f,1,2,100,1e-8)                                          # il secondo e il terzo parametro sono i limiti dell'intervallo di prova (1,2)
scipy_bis_res = optimize.root_scalar(f, method = 'bisect', bracket=[1, 2], maxiter = 100, rtol=1e-8)            #bracket è l'intervallo


print("Soluzione trovata dopo " + str(len(my_bis_res)))                         # Printa le iterazioni fatte dalla mia funzione
print("SOLUTIONS " + str(my_bis_res))                                           # Printa le soluzioni per ogni iterazione
print("MY BISECTION RESULT = " + str(my_bis_res[len(my_bis_res)-1]))            # Printa la soluzione finale trovata
print("SCIPY BISECTION RESULT = ")                                              # Printa la soluzione di SCIPY
print(str(scipy_bis_res))
print("\n")


print("#################################### NEWTON")
my_new_res = newton(f,Df,1,100,1e-8)
scipy_new_res = optimize.root_scalar(f, method = 'newton', x0=1, fprime = Df, maxiter = 100, rtol=1e-8)

print("SOLUTIONS " + str(my_new_res))
print("MY NEWTON RESULT = " + str(my_new_res[len(my_new_res)-1]))
print("SCIPY NEWTON RESULT = ")
print(str(scipy_new_res))
print("\n")


print("#################################### SECANTE")
my_sec_res = secante(f,1,2,100,1e-8)                                            # il secondo e il terzo parametro sono i limiti dell'intervallo di prova (1,2)
scipy_sec_res = optimize.root_scalar(f, x0=1, x1=2, bracket=[1, 2], method = "secant", rtol=1e-8)       #bracket è l'intervallo

print("Soluzione trovata dopo " + str(len(my_sec_res)))
print("SOLUTIONS " + str(my_sec_res))
print("MY SECANT RESULT = " + str(my_sec_res[len(my_sec_res)-1]))
print("SCIPY SECANT RESULT = ")
print(str(scipy_sec_res))
print("\n")


plt.plot(np.arange(len(my_bis_res)), my_bis_res, color='blue', linewidth=3)     # Grafico Bisezione
plt.plot(np.arange(len(my_new_res)), my_new_res, color='green', linewidth=3)    # Grafico Newton
plt.plot(np.arange(len(my_sec_res)), my_sec_res, color='red', linewidth=3)      # Grafico Secante
plt.show()
