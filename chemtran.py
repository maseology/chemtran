
# analytical solutions to 1D solute convective-dispersive transport through porous media with linear adsorption, zero-order production and first-order decay
# ref: van Genuchten, M.Th., 1981. Analytical solutions for chemical transport with simultaneous adsorptions, zero-order production and first-order decay. Journal of Hydrology 49. 213--233.
#
# Assumptions:
# - 1D
# - single-ion solute
# - semi-infinite domain

import math

class chemtran:  
    Ci = 0.0
    Cb = 0.0

    def __init__(self,D,v,R,mu):
        self.D = D # dispersion coefficient
        self.v = v # pore water velocity
        self.R = R # retardation factor
        self.mu = mu # general first-order decay constant

    def muFromk(self,k):
        # k is the first-order decay coefficient; zero when not needed (alpha in van Genuchten, 1981); assumes alpha=beta
        return k*self.R # [1/day] (pg.216 below eq.5)


    ## Solutions for a third-type (flux) boundary condition

    ### Case A1:
    ###  c(x,0) = Ci
    ###  flux = if 0<t<=t0 vC0 = c(0,t) elseif t>t0 0
    ###  \frac{\partial c}{\partial x}(\infty, t)=0
    def cxtA1(self, x, t, C0, t0, g):    
        # C0 concentration at surface boundary [meq/L]
        # t0 duration of concentration pulse [days]
        # g zero-order liquid-phase decay constant (gamme in van Genuchten, 1981) [meq/L/day]
        if self.mu != 0:
            # eq.10
            c = (C0-g/self.mu)*self.__A(x,t,self.v,self.D,self.R,self.mu) + self.__B(x,t,self.Ci,self.v,self.D,self.R,self.mu, g)
            if t>t0: c -= C0*self.__A(x,t-t0,self.v,self.D,self.R,self.mu) 
            return c
        else:
            # eq.c-3
            c = self.Ci + (C0-self.Ci)*self.__U(x,t,self.v,self.D,self.R) + self.__V(x,t,self.v,self.D,self.R, g)
            if t>t0: c -= C0*self.__U(x,t-t0,self.v,self.D,self.R)
            return c

    ### Case A2:
    ###  t, t0 --> infinity (steady-state)
    ###  c(x,0) = Ci
    ###  flux = vC0 = c(0,inf)
    ###  \frac{\partial c}{\partial x}(\infty)=0
    def cxA2(self, x, C0, g): # eq.14
        u = self.__u()
        return g/self.mu + (C0-g/self.mu)*(2*self.v/(self.v+u))*math.exp((self.v-u)*x/2/self.D)

    ### Case A3:
    ###  same as equation 10, but with a background concentration Cb build using a steady-state solution (eq.17)
    def cxiA3(self, x, Cb, g): # eq.17
        u = self.__u()
        return g/self.mu + (Cb-g/self.mu)*(2*self.v/(self.v+u)) * math.exp((self.v-u)*x/2/self.D)

    def cxtA3(self, x, t, Cb, C0, t0, g):
        if self.mu != 0:           
            # eq.18
            c = (C0-Cb)*self.__A(x,t,self.v,self.D,self.R,self.mu) + self.__E(x,Cb,self.v,self.D,self.mu, g)
            if t>t0: c -= Cb*self.__A(x,t-t0,self.v,self.D,self.R,self.mu) # Note: there was a typo in the original paper C0 should have benn Cb
            return c
        else:
            pass

    ### Case A4:
    ###  c(x,0) = Ci
    ###  flux = vC0exp(-lt) # same as A1 except a pulse-type boundary condition is replaced by an exponentially decaying source
    ###  \frac{\partial c}{\partial x}(\infty, t)=0
    def cxtA4(self, x, t, C0, l, g): # eq.20
        if self.mu != l*self.R:
            return C0*self.__F(x,t,self.v,self.D,self.R,self.mu, l) + self.__B(x,t,self.Ci,self.v,self.D,self.R,self.mu, g) - g/self.mu*self.__A(x,t,self.v,self.D,self.R,self.mu)
        else:
            return C0*self.__G(x,t,self.v,self.D,self.R,self.mu) + self.__B(x,t,self.Ci,self.v,self.D,self.R,self.mu, g) - g/self.mu*self.__A(x,t,self.v,self.D,self.R,self.mu)       

    ### Case A5:
    ###  (combination of A3 initial state, and A4 decay soltuion)
    def cxtA5(self, x, t, Cb, C0, l, g): # eq.20
        if self.mu != l*self.R:
            return C0*self.__F(x,t,self.v,self.D,self.R,self.mu, l) - Cb*self.__A(x,t,self.v,self.D,self.R,self.mu) + self.__E(x,Cb,self.v,self.D,self.mu, g)
        else:
            return C0*self.__G(x,t,self.v,self.D,self.R,self.mu) - Cb*self.__A(x,t,self.v,self.D,self.R,self.mu) + self.__E(x,Cb,self.v,self.D,self.mu, g)


    ## Solutions for a third-type (flux) boundary condition

    ### Case B1:
    def cxtB1(self):
        pass



    # extra functions
    def __u(self):
        return self.v*math.sqrt(1+4*self.mu*self.D/self.v**2) # eq.13
    def __w(self, l):
        return self.v*math.sqrt(1+4*self.D*(self.mu-l*self.R)/self.v**2) # eq.22

    # Analytical functions
    def __A(self, x, t, v, D, R, mu): # eq.11
        u = self.__u()
        r =  v/(v+u) * math.exp((v-u)*x/2/D) * math.erfc((R*x-u*t)/2/math.sqrt(D*R*t))
        r += v/(v-u) * math.exp((v+u)*x/2/D) * math.erfc((R*x+u*t)/2/math.sqrt(D*R*t))
        r += v**2/2/mu/D * math.exp(v*x/D-mu*t/R) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t))
        return r
    def __B(self, x, t, Ci, v, D, R, mu, g): # eq.12
        r =  math.erfc((R*x-v*t)/2/math.sqrt(D*R*t))/2 + math.sqrt(v**2*t/math.pi/R/D) * math.exp(-(R*x-v*t)**2/4/D/R/t)
        r -= (1+v*x/D+v**2*t/D/R)/2 * math.exp(v*x/D) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t))
        r *= (g/mu-Ci) * math.exp(-mu*t/R)
        r += g/mu + (Ci-g/mu)*math.exp(-mu*t/R)
        return r    

    def __E(self, x, Cb, v, D, mu, g): # eq.17
        u = self.__u()
        return g/mu + (Cb-g/mu)*(2*v/(v+u)) * math.exp((v-u)*x/2/D)

    def __F(self, x, t, v, D, R, mu, l): # eq.21
        w = self.__w(l)
        r =  v/(v+w) * math.exp((v-w)*x/2/D) * math.erfc((R*x-w*t)/2/math.sqrt(D*R*t))
        r += v/(v-w)/x * math.exp((v+w)*x/2/D) * math.erfc((R*x+w*t)/2/math.sqrt(D*R*t))
        r *= math.exp(-l*t)
        r += v**2/2/D/(mu-l*R) * math.exp(v*x/D-mu*t/R) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t))
        return r
    def __G(self, x, t, v, D, R, mu): # eq.23
        r =  math.erfc((R*x-v*t)/2/math.sqrt(D*R*t)) / 2
        r += math.sqrt(v**2*t/math.pi/D/R) * math.exp(-(R*x-v*t)**2/4/D/R/t)
        r -= (1+v*x/D+v**2*t/D/R) * math.exp(v*x/D) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t)) / 2
        r *= math.exp(-mu*t/R)

    def __U(self, x, t, v, D, R): # eq.c-4
        r =  math.erfc((R*x-v*t)/2/math.sqrt(D*R*t))/2 + math.sqrt(v**2*t/math.pi/D/R) * math.exp(-(R*x-v*t)**2/4/D/R/t)
        r -= (1+v*x/D+v**2*t/D/R)/2 * math.exp(v*x/D) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t))
        return r
    def __V(self, x, t, v, D, R, g): # eq.c-5
        r =  t - (t/2-R*x/2/v-D*R/2/v**2)*math.erfc((R*x-v*t)/2/math.sqrt(D*R*t))
        r -= math.sqrt(t/4/math.pi/D/R) * (R*x+v*t+2*D*R/v) * math.exp(-(R*x-v*t)**2/4/D/R/t)
        r += (t/2-D*R/2/v**2+(R*x+v*t)**2/4/D/R) * math.exp(v*x/D) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t))
        r *= g/R
        return r    