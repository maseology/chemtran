
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
        # g zero-order liquid-phase decay constant (gamma in van Genuchten, 1981) [meq/L/day]
        if self.mu != 0:
            # eq.10
            c = (C0-g/self.mu)*self.__A(x,t,self.v,self.D,self.R,self.mu) + self.__B(x,t,self.Ci,self.v,self.D,self.R,self.mu, g)
            if t>t0: c -= C0*self.__A(x,t-t0,self.v,self.D,self.R,self.mu) 
            return c
        else:
            # eq.c-3 (Case C1)
            c = self.Ci + (C0-self.Ci)*self.__U(x,t,self.v,self.D,self.R) + self.__V(x,t,self.v,self.D,self.R, g)
            if t>t0: c -= C0*self.__U(x,t-t0,self.v,self.D,self.R)
            return c

    ### Case A2:
    ###  t, t0 --> infinity (steady-state)
    ###  c(x,0) = Ci
    ###  flux = vC0 = c(0,inf)
    ###  \frac{\partial c}{\partial x}(\infty)=0
    def cxA2(self, x, C0, g): # eq.14
        if self.mu != 0:
            u = self.__u()
            return g/self.mu + (C0-g/self.mu)*(2*self.v/(self.v+u))*math.exp((self.v-u)*x/2/self.D)
        else:
            # eq.c-6 (Case C2)
            return C0 + g*(self.v*x+self.D)/self.v**2

    ### Case A3:
    ###  same as equation 10, but with a background concentration Cb build using a steady-state solution (eq.17)
    def cxiA3(self, x, g): # eq.17 (initial conditions)
        if self.mu != 0:
            return self.__E(x,self.Cb,self.v,self.D,self.mu, g)
        else:
            # eq.c-7 (Case C3)
            return self.Cb + g*(self.v*x+self.D)/self.v**2
            


    def cxtA3(self, x, t, C0, t0, g):
        if self.mu != 0:           
            # eq.18
            c = (C0-self.Cb)*self.__A(x,t,self.v,self.D,self.R,self.mu) + self.__E(x,self.Cb,self.v,self.D,self.mu, g)
            if t>t0: c -= self.Cb*self.__A(x,t-t0,self.v,self.D,self.R,self.mu) # Note: there was a typo in the original paper C0 should have been Cb
            return c
        else:
            # eq.c-8 (Case C3)
            c = self.Cb + (C0-self.Cb)*self.__U(x,t,self.v,self.D,self.R) + g*(self.v*x+self.D)/self.v**2
            if t>t0: c -= self.Cb*self.__U(x,t-t0,self.v,self.D,self.R) # Note: (assuming) there was a typo in the original paper C0 should have been Cb (as found in A3)
            return c

    ### Case A4:
    ###  c(x,0) = Ci
    ###  flux = vC0exp(-lt) # same as A1 except a pulse-type boundary condition is replaced by an exponentially decaying source
    ###  \frac{\partial c}{\partial x}(\infty, t)=0
    def cxtA4(self, x, t, C0, l, g): # eq.20
        if self.mu != 0:
            if self.mu != l*self.R:
                return C0*self.__F(x,t,self.v,self.D,self.R,self.mu, l) + self.__B(x,t,self.Ci,self.v,self.D,self.R,self.mu, g) - g/self.mu*self.__A(x,t,self.v,self.D,self.R,self.mu)
            else:
                return C0*self.__G(x,t,self.v,self.D,self.R,self.mu) + self.__B(x,t,self.Ci,self.v,self.D,self.R,self.mu, g) - g/self.mu*self.__A(x,t,self.v,self.D,self.R,self.mu)       
        else:
            # eq.c-9 (Case C4)
            return self.Ci - self.Ci*self.__U(x,t,self.v,self.D,self.R) + C0*self.__W(x,t,self.v,self.D,self.R,l) + self.__V(x,t,self.v,self.D,self.R, g)

    ### Case A5:
    ###  (combination of A3 initial state, and A4 decay soltuion)
    def cxtA5(self, x, t, C0, l, g): # eq.20
        if self.mu != 0:
            if self.mu != l*self.R:
                return C0*self.__F(x,t,self.v,self.D,self.R,self.mu, l) - self.Cb*self.__A(x,t,self.v,self.D,self.R,self.mu) + self.__E(x,self.Cb,self.v,self.D,self.mu, g)
            else:
                return C0*self.__G(x,t,self.v,self.D,self.R,self.mu) - self.Cb*self.__A(x,t,self.v,self.D,self.R,self.mu) + self.__E(x,self.Cb,self.v,self.D,self.mu, g)
        else:
            # eq.c-12 (Case C5)
            return self.Cb - self.Cb*self.__U(x,t,self.v,self.D,self.R) + C0*self.__W(x,t,self.v,self.D,self.R,l) + g*(self.v*x+self.D/self.v**2)

    ## Solutions for a third-type (flux) boundary condition

    ### Case B1:
    ###  c(x,0) = Ci
    ###  c(0,t) = if 0<t<=t0 C0 elseif t>t0 0
    ###  \frac{\partial c}{\partial x}(\infty, t)=0
    def cxtB1(self, x, t, C0, t0, g): 
        if self.mu != 0:   
            # C0 concentration at surface boundary [meq/L]
            # t0 duration of concentration pulse [days]
            # g zero-order liquid-phase decay constant (gamma in van Genuchten, 1981) [meq/L/day]
            # eq. 26
            c = (C0-g/self.mu)*self.__H(x,t,self.v,self.D,self.R) + self.__M(x,t,self.Ci,self.v,self.D,self.R,self.mu,g) #, Ci, v, D, R, mu, g
            if t>t0: c -= C0*self.__H(x,t-t0,self.v,self.D,self.R)
            return c
        else:
            # equation c-13 (Case D1)
            c = (self.Ci + (C0-self.Ci)*self.__X(x,t,self.v,self.D,self.R) + self.__Y(x,t,self.v,self.D,self.R,g))
            if t>t0: c -= C0*self.__X(x,t-t0,self.v,self.D,self.R)
            return c

    ### Case B2:
    ###  t, t0 --> infinity (steady-state)
    ###  c(0) = C0
    ###  \frac{\partial c}{\partial x}(\infty)=0    
    def cxB2(self, x, C0, g): # eq.29
        if self.mu != 0:
            u = self.__u()
            return g/self.mu + (C0-g/self.mu)*math.exp((self.v-u)*x/2/self.D)   
        else:
            return C0 + g*x/self.v # eq.c-16 (Case D2)

    ### Case B3:
    ###  similar to equation A3, but with a background concentration Cb build using a steady-state solution (eq.30)
    def cxiB3(self, x, g): # eq.30 (initial conditions)
        if self.mu != 0:
            return self.__N(x,self.Cb,self.v,self.D,self.mu,g)
        else:
            return self.Cb + g*x/self.v # eq.c-17 (Case D3)

    def cxtB3(self, x, t, C0, t0, g):
        if self.mu != 0:           
            # eq.18
            c = (C0-self.Cb)*self.__H(x,t,self.v,self.D,self.R) + self.__N(x,self.Cb,self.v,self.D,self.mu,g)
            if t>t0: c -= self.Cb*self.__H(x,t-t0,self.v,self.D,self.R) # Note: (assuming) there was a typo in the original paper C0 should have been Cb (as found in A3)
            return c
        else:
            # eq.c-18 (Case D3)
            c = self.Cb + (C0-self.Cb)*self.__X(x,t,self.v,self.D,self.R) + g*x/self.v
            if t>t0: c -= C0*self.__X(x,t-t0,self.v,self.D,self.R)


    ### Case B4:
    ###  same as A4, only with exponentially decaying boundary condition: c(0,t) = C0*exp(-lambda*t)
    ###  c(x,0) = Ci
    ###  \frac{\partial c}{\partial x}(\infty, t)=0
    def cxtB4(self, x, t, C0, l, g): # eq.33
        if self.mu != 0:
            return C0*self.__P(x,t,self.v,self.D,self.R,l) + self.__M(x,t,self.Ci,self.v,self.D,self.R,self.mu,g) - g/self.mu*self.__H(x,t,self.v,self.D,self.R)
        else:
            # eq.c-19 (Case D4)
            return self.Ci - self.Ci*self.__X(x,t,self.v,self.D,self.R) + C0*self.__Z(x,t,self.v,self.D,self.R,l) + self.__Y(x,t,self.v,self.D,self.R,g)

    ### Case B5:
    ###  analogous to A5, steady-state type initinal concentration (as in B3), an exponentially decaying boundary condition and for a semi-infintite medium 
    def cxtB5(self, x, t, C0, l, g): # eq.35
        if self.mu != 0:
            return C0*self.__P(x,t,self.v,self.D,self.R,l) - self.Cb*self.__H(x,t,self.v,self.D,self.R) + self.__N(x,self.Cb,self.v,self.D,self.mu,g)
        else:
            # eq.c-21(Case D5)
            return self.Cb - self.Cb*self.__X(x,t,self.v,self.D,self.R) + C0*self.__Z(x,t,self.v,self.D,self.R,l) + g*x/self.v



    ##################################
    ####  supplemental functions  ####
    ##################################
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
        return r

    def __H(self, x, t, v, D, R): # eq.27
        u = self.__u()
        r =  math.exp((v-u)*x/2/D) * math.erfc((R*x-u*t)/2/math.sqrt(D*R*t)) / 2
        r += math.exp((v+u)*x/2/D) * math.erfc((R*x+u*t)/2/math.sqrt(D*R*t)) / 2
        return r
    def __M(self, x, t, Ci, v, D, R, mu, g) : # eq.28
        r =  math.erfc((R*x-v*t)/2/math.sqrt(D*R*t)) / 2
        r += math.exp(v*x/D) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t)) / 2
        r *= (g/mu-Ci) * math.exp(-mu*t/R)
        r += g/mu + (Ci-g/mu) * math.exp(-mu*t/R)
        return r

    def __N(self, x, Cb, v, D, mu, g): # eq.30
        u = self.__u()
        r =  g/mu
        r += (Cb-g/mu) * math.exp((v-u)*x/2/D)

    def __P(self, x, t, v, D, R, l):
        w = self.__w(l)
        r =  math.exp((v-w)*x/2/D) * math.erfc((R*x-w*t)/2/math.sqrt(D*R*t)) / 2
        r += math.exp((v+w)*x/2/D) * math.erfc((R*x+w*t)/2/math.sqrt(D*R*t)) / 2
        r *= math.exp(-l*t)
        return r

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

    def __W(self, x, t, v, D, R, l):
        s = v*math.sqrt(1-4*l*D*R/v**2) # eq.c-11
        r =  math.exp((v-s)*x/2/D) * math.erfc((R*x-s*t)/2/math.sqrt(D*R*t)) * v/(v+s)
        r += math.exp((v+s)*x/2/D) * math.erfc((R*x+s*t)/2/math.sqrt(D*R*t)) * v/(v-s)
        r *= math.exp(-l*t)
        r -= math.exp(v*x/D) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t)) * v**2/2/l/D/R
        return r

    def __X(self, x, t, v, D, R): # eq.c-14
        r =  math.erfc((R*x-v*t)/2/math.sqrt(D*R*t)) / 2
        r += math.exp(v*x/D) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t)) / 2
        return r

    def __Y(self, x, t, v, D, R, g): # eq.c-15
        r =  t + (R*x-v*t)/2/v * math.erfc((R*x-v*t)/2/math.sqrt(D*R*t))
        r -= (R*x+v*t)/2/v * math.exp(v*x/D) * math.erfc((R*x+v*t)/2/math.sqrt(D*R*t))
        r *= g/R
        return r

    def __Z(self, x, t, v, D, R, l): # eq.c-20
        s = v*math.sqrt(1-4*l*D*R/v**2) # eq.c-11
        r =  math.exp((v-s)*x/2/D) * math.erfc((R*x-s*t)/2/math.sqrt(D*R*t)) / 2
        r += math.exp((v+s)*x/2/D) * math.erfc((R*x+s*t)/2/math.sqrt(D*R*t)) / 2
        r *= math.exp(-l*t)
        return r