'''Adapted from Chapter 1.2 The Nonlinear Pendulum of
    Kinzel and Reents.
    Coded by: Christian Alis
    Annotated by: May Lim
    
    Shows Figs 1.1 - 1.4
'''

from scipy.special import ellipk, ellipj
from sympy.parsing.sympy_parser import parse_expr
from numpy.fft import fft
import cProfile
import numpy as np
import matplotlib.pyplot as plt
import sympy

print "\nNonlinear Pendulum\n"

#Eq 1.7 for period T, with \sqrt(l/g) = 1
def T(phi0):
    return 4*ellipk(np.sin(phi0/2)**2)
#Read help(ellipk) -- crucial info for phi0 = np.pi
#... and for figuring out the convention that is in use

#Alternative expression for Eq 1.7
#T = lambda phi0: 4*ellipk(np.sin(phi0/2)**2)

#Generate data for Fig 1.1 - period T as a function of amplitude phi0
phis = np.linspace(0, 3.14, 200)
ts = T(phis)

#Plotting routines for Fig 1.1: take note of...
# ... launch a new figure
plt.figure(figsize=(5,4))
plt.plot(phis/np.pi, ts/np.pi)
plt.ylim(0, 30/np.pi)
plt.xlabel(r'$\phi_0$', size='xx-large')
plt.ylabel('T', rotation="vertical", size='xx-large')
#... tick labels
plt.xticks([0, 0.25, 0.5, 0.75, 1], ["0", "$\pi$/4", "$\pi$/2", "3$\pi$/4", "$\pi$"])
plt.yticks([0, 2, 4, 6, 8, 10], ["0", "2$\pi$", "4$\pi$", "6$\pi$", "8$\pi$", "10$\pi$"])
plt.savefig('Fig(1.1).pdf')
plt.show()


####QUESTIONS####
##How should the code be changed if I want the Greek symbol \pi instead?
'''I used "$\pi$" instead of "Pi". '''

##What other plot options for plt.xlabel are available?
'''The plt.xlabel set the x axis label. This function offers other option like 
changing the fontsize, fontstyle and vertical and horizontal alignment'''

#What would happen if plt.show() is not commented out? Try it.
'''Figure 1.1 is shown first before the equation.'''

#The function t(phi) is obtained from Eq 1.6, but we want to calculate
#the inverse phi(t) -- this is also available
#Check out: np.info(ellipj)
def sinuspsi(t, phi0):
    return ellipj(t, np.sin(phi0/2)**2)[0] 
#Alternatively,
#sinuspsi = lambda t, phi0: ellipj(t, np.sin(phi0/2)**2)[0]
#More info on the Jacobi elliptic function: 
#http://keisan.casio.com/has10/SpecExec.cgi?id=system/2006/1180573437
#Note the first element in the output of ellipj above is sn(u,k) is
#the left-hand term of Eq 1.5

#Eq 1.5 retrieves phi. Normalize by dividing with phi0
def phinorm(x, phi0):
    return 2*np.arcsin(np.sin(phi0/2) * sinuspsi(x*T(phi0), phi0))/phi0
#Alternatively,
#phinorm = lambda x, phi0: 2*np.arcsin(np.sin(phi0/2) * sinuspsi(x*T(phi0), phi0))/phi0

#Data/Plotting routines for Fig 1.2
phi0s = np.array([0.1, 0.8, 0.95, 0.99, 0.999]) * np.pi
tnorms = np.linspace(0,1)
phinorms = [phinorm(tnorms, phi0) for phi0 in phi0s]
plt.figure(figsize=(5,4))
plt.plot(tnorms.T, np.transpose(phinorms))
plt.axhline(0)
plt.xlabel('t/T')
plt.ylabel('phi/phi0')
plt.savefig('Fig(1.2).pdf')

####ACTIVITY/QUESTIONS####
#Double-check the expressions sinuspsi and phinorm
#Where is the zero reference point of phi?
'''The zero reference point of phi is the stable of equilibrium of the pendulum. '''

#When phi/phi0 = 1, where is the bob?
'''The bob is at the maximum angle. '''

#On the programming side, you may want to read about list comprehension 
#here http://www.python.org/dev/peps/pep-0202/

#Describe the anharmonic oscillations for an amplitude phi0 = 0.999*pi
'''The oscillations become square wave.'''

#Data and plotting routines for Fig 1.3
list_ = phinorm(np.arange(0, 1, 0.01), 0.999*np.pi)
foulist = (np.abs(fft(list_))/np.sqrt(len(list_)))[:15]
plt.figure(figsize=(5,4))
plt.scatter(range(15), foulist)
plt.axhline(0)
plt.xlim(0, 15)
plt.ylim(-1, 8)
plt.xlabel("s")
plt.ylabel("Abs[b(s)]")
plt.savefig('Fig(1.3).pdf')

#Fig 1.4: phase-space plot of phidot vs. phi
#This code mirrors the symbolic source code of Mathematica
m = sympy.Symbol("m")
psi = sympy.Symbol("psi")
phi0 = sympy.Symbol("phi0")
f = 1/sympy.sqrt(1 - m*sympy.sin(psi)**2)
g = f.series(m, 0, 10)
tseries = 4*sympy.integrate(g, (psi, 0, sympy.pi/2))
print tseries.subs(m, sympy.sin(phi0/2)**2)

#phi0 = sympy.Symbol('%f' %phi0s[0])
#m = sympy.sin(phi0/2)**2
#f = 1/sympy.sqrt(1 - m*sympy.sin(psi)**2)
#g = f.series(m, 0, 10)
#tseries = 4*sympy.integrate(g, (psi, 0, sympy.pi/2))
#tseriesval = parse_expr(str(tseries)[0:194])
#print float(tseriesval)
#tseries1 = 4.0*ellipk((np.sin(float(phi0s[0])/2.0))**2)
#print tseries1
####ACTIVITY####
#Rewrite the code for Fig 1.4 so that the integration will be done numerically
'''Doing the integral analytically'''
#let phi0 = np.pi/6.0
def an(phi0):
    return ellipk(np.sin(phi0/2)**2)
print 'Analitically, tseries = ' + str(an(np.pi/6.0))
'''Doing the integral numerically using trapezoidal rule'''
def trap(phi0):
    m = phi0/2.0
    def f(psi):
        return 1.0/np.sqrt(1 - m*np.sin(psi)**2)
    
    N = 10
    a = 0.0
    b = np.pi/2.0
    h = (b-a)/N
    
    s = 0.5*f(a) + 0.5*f(b)
    for k in range(1,N):
        s += f(a+k*h)
    return h*s
print 'Computationally, tseries = ' + str(trap(np.pi/6.0))

#Profile both versions - see http://docs.python.org/library/profile.html
cProfile.run('an(np.pi/6.0)')
cProfile.run('trap(np.pi/6.0)')
#Make a contour plot
def e(phidot, phi):
    return phidot**2/2. - np.cos(phi)
#Alternatively,
#e = lambda phidot, phi: phidot**2/2. - np.cos(phi)
#mgrid is an optimal way of creating a grid for evaluating a function
y, x = np.mgrid[-3:3:0.1, -np.pi:np.pi:0.1]
z = e(y, x)
plt.figure(figsize=(5,4))
plt.contour(x/np.pi, y, z, levels=[-0.5, 0, 0.5, 1, 1.5, 2])
plt.xlabel("phi")
plt.ylabel("phidot")
plt.xticks([-1, -0.5, 0, 0.5, 1], ["-$\pi$", "-$\pi$/2", 0, "$\pi$/2", "$\pi$"])
plt.savefig('Fig(1.4).pdf')
plt.show()


####EXTRAS####
#Save the figures to .pdf files (automatically)
#Save the figures with specific dimensions (e.g. width = 3 inches)