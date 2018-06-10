from __future__ import print_function, division
from six.moves import range

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

import sys

"""
Two equations system

d(rhoe)/dt + d(rhoe*vh)/dx = d(D/Te*d(pe)/dx)/dx
d(rhoe*ee)/dt + d(rhoe*ee*vh)/dx = -pe*d(vh)/dx + d(L*d(Te)/dx + gamma/(gamma-1)*D*d(pe)/dx)/dx

with p = (gamma-1) rho e
"""

def printProgress (iteration, total, prefix = '', suffix = '', decimals = 1,  barLength = 100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    formatStr       = '{0:' + str(decimals+4) + '.' + str(decimals) + 'f}'
    percents        = formatStr.format(100 * (iteration / float(total)))
    filledLength    = int(round(barLength * iteration / float(total)))
    bar             = '-' * filledLength + ' ' * (barLength - filledLength)
    print('\r{0:s} |{1:s}| {2:s}% {3:s}'.format(prefix, bar, percents,  suffix), end='', file=sys.stdout, flush=True)
    if iteration == total:
        print('', end = '\n', file=sys.stdout, flush=True)

def flux_num_LF(S):
    """
    numerical flux with a Lax-Fridrichs scheme for the convective part
    """
    ## F_{j+1/2} = 0.5*( F(U_{j})+F(U_{j+1})+(U_j-U_{j+1})*dx/dt )
    S.F_conv[:] =  0.5*( S.flux_conv(S.t,S.X_glob[:-1],S.sol_glob[:-1,:]) + \
                        S.flux_conv(S.t,S.X_glob[1:],S.sol_glob[1:,:]) + \
                        (S.sol_glob[:-1,:]-S.sol_glob[1:,:])*S.dx/S.dt )

def flux_num_UW(S):
    """
    numerical flux with a upwind scheme for the convective part
    """
    S.F_conv[:] =  S.flux_conv(S.t,S.X_glob[:-1],S.sol_glob[:-1,:])

def nc_num_MW(S):
    """
    numerical contribution of the non conservative term
    with the scheme proposed by Massot and Wargnier
    """
    dummy = S.P.gamma - 1.
    Terond = dummy * S.sol[:,1] / S.sol[:,0]
    pei = dummy * .5 * ( S.sol_glob[1:,1] + S.sol_glob[:-1,1] )
    Tei = dummy * .5 * ( S.sol_glob[1:,1]/S.sol_glob[1:,0] + S.sol_glob[:-1,1]/S.sol_glob[:-1,0] )
    gpei = dummy / S.dx * (S.sol_glob[1:,1] - S.sol_glob[:-1,1])
    termegradpe = (Tei[1:] - Terond) / Tei[1:] * gpei[1:] \
        - (Tei[:-1] - Terond) / Tei[:-1] * gpei[:-1]
    v = S.exsol.vh(S.t, S.X_int)
    rhoev = S.F_conv[:, 0]
    rhoeeev = S.F_conv[:, 1]
    dummy = S.P.gamma / (S.P.gamma-1.)
    v = S.exsol.vh(S.t, S.X)
    
    grhoeee = np.ediff1d(S.sol_glob[:-1,1])
    rhoeeeR, rhoeeeL = S.P.peR/(S.P.gamma-1), S.P.peL/(S.P.gamma-1)
    if self.P.scheme_conv == "UW":
        seuil = .35*abs(rhoeeeR-rhoeeeL)
    elif self.P.scheme_conv == "LF":
        seuil = .17*abs(rhoeeeR-rhoeeeL)
    else:
        seuil = 0
    grhoeee[np.where(grhoeee>seuil)] = seuil
    grhoeee[np.where(grhoeee<-seuil)] = -seuil
    gTe = np.ediff1d(S.sol_glob[:-1,1] / S.sol_glob[:-1,0])
    TeR, TeL = S.P.peR/S.P.rhoeR, S.P.peL/S.P.rhoeL
    seuil = .1*abs(TeR-TeL)
    gTe[np.where(gTe>seuil)] = seuil
    gTe[np.where(gTe<-seuil)] = -seuil
    
    correction = (S.P.gamma-1.) * v * grhoeee - S.P.gamma * v * S.sol[:,0] * gTe
    S.nc[:,1] = dummy * Terond * np.ediff1d(rhoev) - np.ediff1d(rhoeeev) \
        - correction #+ dummy * S.P.diffD * termegradpe

def nc_num_df(S):
    """
    numerical contribution of the non conservative term
    with a finite differences approximation
    """
    divv = np.ediff1d(S.exsol.vh(S.t, S.X_int))
    S.nc[:, 1] = (S.P.gamma-1.) * S.sol[:, 1] * divv

def diff_num_centered(S):
    """
    numerical flux with a centered scheme for the dissipative part
    """
    # Ui = [rhoe, rhoe ee] at the interfaces
    Ui = .5*(S.sol_glob[1:,:] + S.sol_glob[:-1,:])
    # dUidx = [grad Te, grad pe] at the interfaces
    dUidx = np.stack((
        (S.P.gamma-1.)/S.dx *  np.diff(S.sol_glob[:,1] / S.sol_glob[:,0]),
        (S.P.gamma-1.)/S.dx *  np.diff(S.sol_glob[:,1])
        ), axis = -1)
    S.F_diff[:] = S.flux_diff(S.t + S.t_diff, S.X, Ui, dUidx)

class exacte_sol():
    """
    Exact solution for the reduced system on the electronic variables
    Case of a Riemann problem initialization
    The progressive waves corresponding to the shocks are computed through the exact formulae
    """
    def __init__(self, P):
        self.P = P
        if (self.P.onde==1):
            self._init_onde1()
        elif (self.P.onde==3):
            self._init_onde3()
        else:
            print("Error: choose 1 or 3 to initialize wave")
        if self.P.sigma > 0:
            self.xdiscont = P.Lx*0.25
        elif self.P.sigma < 0:
            self.xdiscont = P.Lx*0.75
        else:
            self.xdiscont = P.Lx*0.5

    def _init_onde3(self):
        # compute the left state from the right state and the left total pressure
        self.P.mh =  -np.sqrt( self.P.rhohR*(self.P.pL+self.P.mu2*self.P.pR)/(1.-self.P.mu2) )
        self.P.vhL = self.P.vhR - (self.P.pL-self.P.pR)/self.P.mh
        self.P.sigma = self.P.vhR - self.P.mh/self.P.rhohR
        self.P.uhL = self.P.vhL - self.P.sigma
        self.P.uhR = self.P.vhR - self.P.sigma
        self.P.rhohL = self.P.mh/self.P.uhL
        self.P.rhoeL = self.P.rhoeR*self.P.uhR/self.P.uhL
        self.P.peL = self.P.peR*self.P.uhR/( self.P.gamma*self.P.uhL-(self.P.gamma-1)*self.P.uhR )
        # compute the regularization
        d = self.P.diffD/(self.P.vhR-self.P.sigma)
        l = self.P.diffL/(self.P.vhR-self.P.sigma)
        dummy = 1.0/d+self.P.gamma/(self.P.gamma-1.)*self.P.rhoeR/l 
        dd    = dummy**2 - 4.0/(self.P.gamma-1.0)*self.P.rhoeR/(l*d)
        lambdap = 0.5*(dummy+np.sqrt(dd)) # first eigenvalue
        lambdam = 0.5*(dummy-np.sqrt(dd)) # second eigenvalue
        self.lambdam = lambdam
        self.lambdap = lambdap
        self.d = d
        self.P.TeR = self.P.peR / self.P.rhoeR
        R = (1./(self.P.vhL-self.P.vhR) - self.P.gamma/(self.P.sigma-self.P.vhR))/(self.P.gamma-1.)
        self.ypo = - (1.-self.P.gamma*d*lambdam)*self.P.TeR/(R*self.P.diffD*(lambdap-lambdam)*(self.P.gamma-1.))
        self.ymo =   (1.-self.P.gamma*d*lambdap)*self.P.TeR/(R*self.P.diffD*(lambdap-lambdam)*(self.P.gamma-1.))

    def _init_onde1(self):
        # compute the right state
        self.P.mh =  0.
        self.P.vhR = 0.
        self.P.sigma = 0.
        self.P.uhR = 0.
        self.P.uhL = 0.
        self.P.rhohR = 0.
        self.P.rhoeR = 0.
        self.P.peR = 0.
        print("1-wave not yet implemented !!!")

    def vh(self, t, x):
        # the heavy particles velocity is a two states function (discontinuity)
        return self.P.vhL+(self.P.vhR-self.P.vhL)*(x-self.P.sigma*t>self.xdiscont)

    def rhoe(self, t, x):
        xx = x - self.P.sigma*t - self.xdiscont
        Xtr = xx * (xx>0)
        pe = self.P.peR + self.P.rhoeR * (self.ypo*np.exp(self.lambdap*Xtr) + self.ymo*np.exp(self.lambdam*Xtr))
        Te = pe/self.P.rhoeR -self.d*self.lambdap*self.ypo*np.exp(self.lambdap*Xtr) -self.d*self.lambdam*self.ymo*np.exp(self.lambdam*Xtr)
        return pe/Te

    def pe(self, t, x):
        xx = x - self.P.sigma*t - self.xdiscont
        Xtr   = xx * (xx>0)
        pe = self.P.peR + self.P.rhoeR * (self.ypo*np.exp(self.lambdap*Xtr) + self.ymo*np.exp(self.lambdam*Xtr))
        return pe

    def Te(self, t, x):
        xx = x - self.P.sigma*t - self.xdiscont
        Xtr   = xx * (xx>0)
        pe = self.P.peR + self.P.rhoeR * (self.ypo*np.exp(self.lambdap*Xtr) + self.ymo*np.exp(self.lambdam*Xtr))
        Te = pe/self.P.rhoeR -self.d*self.lambdap*self.ypo*np.exp(self.lambdap*Xtr) -self.d*self.lambdam*self.ymo*np.exp(self.lambdam*Xtr)
        return Te

    def rhoeee(self, t, x):
        xx = x - self.P.sigma*t - self.xdiscont
        Xtr   = xx * (xx>0)
        pe = self.P.peR + self.P.rhoeR * (self.ypo*np.exp(self.lambdap*Xtr) + self.ymo*np.exp(self.lambdam*Xtr))
        return 1./(self.P.gamma-1.)*pe


class twoeq():
    def __init__(self, P):
        self.P = P
        if self.P.scheme_conv == "UW":
            self.scheme_conv = flux_num_UW
        else:
            self.scheme_conv = flux_num_LF
        self.scheme_diff = diff_num_centered
        if self.P.scheme_nc == "MW":
            self.scheme_nc = nc_num_MW
        else:
            self.scheme_nc = nc_num_df
        self.P.mu2 = (self.P.gamma-1)/(self.P.gamma+1)
        # mesh
        self.X_glob = np.zeros((self.P.N+2,))
        self.X_glob[1:] = np.linspace(0.0, self.P.Lx, self.P.N+1)
        self.dx = self.X_glob[2]-self.X_glob[1]
        self.X_glob[0] = -self.dx
        self.X_glob += self.dx*0.5
        self.X = self.X_glob[1:-1]
        self.X_int = .5*(self.X_glob[1:] + self.X_glob[:-1])
        # solution
        self.sol_glob = np.zeros((self.P.N+2, 2))
        self.sol = self.sol_glob[1:-1,:]
        # convective fluxes
        self.F_conv = np.zeros((self.P.N+1,2))
        # diffusive fluxes
        self.F_diff = np.zeros((self.P.N+1,2))
        # non conservative term
        self.nc = np.zeros((self.P.N, 2))

        self.t = self.P.t0
        self.exsol = exacte_sol(self.P)
        self.init()

        if self.P.splitting:
            self._onetimestep = self._onetimestep_with_splitting
        else:
            self._onetimestep = self._onetimestep_without_splitting

        if self.P.liste_plot is not None:
            self.liste_plot = self.P.liste_plot
        else:
            self.liste_plot = ['rhoe', 'pe', 'Te']


    def var_rhoe(self):
        return self.sol[:, 0]

    def var_Ee(self):
        return self.sol[:, 1]

    def var_pe(self):
        return (self.P.gamma-1.) * self.sol[:, 1]

    def var_Te(self):
        return (self.P.gamma-1.) * self.sol[:, 1] / self.sol[:, 0]

    def init(self):
        regul = True
        if regul:
            # with regularization
            self.sol[:, 0] = self.exsol.rhoe(0., self.X)      # rhoe
            self.sol[:, 1] = self.exsol.rhoeee(0., self.X)    # rhoe ee
        else:
            # without regulartization
            self.sol[:,0] = self.P.rhoeL + (self.P.rhoeR - self.P.rhoeL)*(self.X>self.xdiscont)                       # rhoe
            self.sol[:,1] = ( self.P.peL + (self.P.peR - self.P.peL)*(self.X>self.xdiscont) )/(self.P.gamma-1)        # rhoe ee

    def _dt_compute_dt(self):
        if self.P.diff_substeps:
            self.dt = min(self.P.tf-self.t,self.P.CFL*self.dx/max(self.P.vhL,self.P.vhR))
        else:
            beta = max(self.P.diffD, (self.P.gamma-1.)/self.P.gamma*self.P.diffL/self.P.rhoeR)
            self.dt = self.P.CFL * .5 / beta * self.dx**2

    def flux_conv(self,t,x,U):
        return np.stack((
            U[:,0]*self.exsol.vh(t,x),
            U[:,1]*self.exsol.vh(t,x)
            ), axis = -1)

    def flux_diff(self, t, x, Ui, dUidx):
        return np.stack((
            self.P.diffD / (self.P.gamma-1.) * Ui[:,0] / Ui[:,1] * dUidx[:,1],
            self.P.diffL * dUidx[:,0] + self.P.gamma * self.P.diffD / (self.P.gamma-1.) * dUidx[:,1]
            ), axis = -1)

    def _boundary(self):
        # Neumann boundary conditions
        self.sol_glob[0,:] = self.sol_glob[1,:]
        self.sol_glob[-1,:] = self.sol_glob[-2,:]

    def _onetimestep_conv(self, split_frac = 1):
        # convective fluxes
        self.scheme_conv(self)
        self.sol -= split_frac * self.dt/self.dx * np.diff(self.F_conv, axis = 0)

    def _onetimestep_diff(self, split_frac = 1):
        # diffusive fluxes
        self.t_diff, nb_iter_diff = 0, 0
        if self.P.diff_substeps:
            # use subteps for diffusion to increase global time step (given by a convective CFL)
            while self.t_diff < split_frac * self.dt:
                nb_iter_diff += 1
                dt_diff = min(split_frac * self.dt - self.t_diff, .5*self.dx**2)
                self.scheme_diff(self)
                self.sol += dt_diff/self.dx * np.diff(self.F_diff, axis = 0)
                self.t_diff += dt_diff
        else:
            #
            self.scheme_diff(self)
            self.sol += split_frac * self.dt/self.dx * np.diff(self.F_diff, axis = 0)

    def _onetimestep_nc(self, split_frac = 1):
        self.scheme_nc(self)
        self.sol -= split_frac * self.dt/self.dx * self.nc

    def _onetimestep_with_splitting(self):
        self._dt_compute_dt()
        if self.dt > 0:
            self._boundary()
            # operators splitting
            self._onetimestep_diff(split_frac = .5)
            self._onetimestep_nc(split_frac = .5)
            self._onetimestep_conv(split_frac = 1)
            self._onetimestep_nc(split_frac = .5)
            self._onetimestep_diff(split_frac = .5)
            self.t += self.dt

    def _onetimestep_without_splitting(self):
        self._dt_compute_dt()
        if self.dt > 0:
            self._boundary()
            # all operators are taken at the same time (no splitting)
            self.t_diff = 0.
            self.scheme_conv(self)
            self.scheme_diff(self)
            self.scheme_nc(self)
            self.sol -=  self.dt/self.dx * np.diff(self.F_conv, axis = 0) \
                - self.dt/self.dx * np.diff(self.F_diff, axis = 0) \
                + self.dt/self.dx * self.nc
            self.t += self.dt


    def run(self):
        print("*"*80)
        print("Begin the simulation")
        dico = {
            'rhoe': (self.var_rhoe, self.exsol.rhoe, r'$\rho_e$'),
            'pe': (self.var_pe, self.exsol.pe, r'$p_e$'),
            'Te': (self.var_Te, self.exsol.Te, r'$T_e$'),
        }
        liste_num = [dico[key][0] for key in self.liste_plot]
        liste_exa = [dico[key][1] for key in self.liste_plot]
        titre = tuple([dico[key][2] for key in self.liste_plot])
        #ax = [plt.subplot(1,len(self.liste_plot), k) for k in range(1, len(self.liste_plot)+1)]
        fig, ax = plt.subplots(1, len(self.liste_plot), figsize = (4*len(self.liste_plot), 4))
        self.lines = []
        for k in range(len(liste_num)):
            #self.lines.append(ax[k].plot(self.X, liste_num[k](), linewidth=1, color='orange', alpha=1)[0])
            self.lines.append(ax[k].scatter(self.X, liste_num[k](), s=2, color='orange', alpha=1))
            self.lines.append(ax[k].plot(self.X, liste_exa[k](self.t, self.X), linewidth=1, color='navy', alpha=1)[0])
            ax[k].set_title(titre[k] + r' at $t={0:4.2f}$'.format(self.t))
            yL, yR = liste_exa[k](0., self.X[0]), liste_exa[k](0., self.X[-1])
            ymin, ymax = min(yL, yR), max(yL, yR)
            dy = .25*(ymax - ymin)
            ax[k].set_ylim([ymin-dy, ymax+dy])
        
        plt.savefig(self.P.dir + "sol_t={:12.10f}.pdf".format(self.t))

        def animate(i):
            ti = min(i*self.P.peraff, self.P.tf)
            while self.t < ti:
                self._onetimestep()
                printProgress(self.t, self.P.tf, prefix = 'Progress:', suffix = 'Complete', barLength =  50)
            for k in range(len(liste_num)):
                #self.lines[2*k].set_data(self.X, liste_num[k]())
                self.lines[2*k].set_offsets(np.c_[self.X, liste_num[k]()])
                self.lines[2*k+1].set_data(self.X, liste_exa[k](self.t, self.X))
                ax[k].set_title(titre[k] + r' at $t={0:4.2f}$'.format(self.t))
            plt.savefig(self.P.dir + "sol_t={:12.10f}.pdf".format(self.t))

        anim = animation.FuncAnimation(fig, animate, frames=100, interval=10)
        plt.show()
        print("\nEnd")
        print("*"*80)
        
    def __str__(self):
        str  = "*"*80 + "\n"
        str += "Informations:\n"
        str += "-"*50 + "\n"
        str += "\tLeft            Right" + "\n"
        str += "-"*50 + "\n"
        str += "p\t{0:10.3e}\t{1:10.3e}\n".format(self.P.pL, self.P.pR)
        str += "rhoh\t{0:10.3e}\t{1:10.3e}\n".format(self.P.rhohL, self.P.rhohR)
        str += "vh\t{0:10.3e}\t{1:10.3e}\n".format(self.P.vhL, self.P.vhR)
        str += "-"*50 + "\n"
        str += "rhoe\t{0:10.3e}\t{1:10.3e}\n".format(self.P.rhoeL, self.P.rhoeR)
        # str += "    \t{0:10.3e}\t{1:10.3e}\n".format(self.sol[0,0], self.sol[-1,0])
        str += "pe\t{0:10.3e}\t{1:10.3e}\n".format(self.P.peL, self.P.peR)
        # str += "    \t{0:10.3e}\t{1:10.3e}\n".format((self.P.gamma-1.) * self.sol[0,1], (self.P.gamma-1.) * self.sol[-1,1])
        str += "Te\t{0:10.3e}\t{1:10.3e}\n".format(self.P.peL/self.P.rhoeL, self.P.peR/self.P.rhoeR)
        # str += "    \t{0:10.3e}\t{1:10.3e}\n".format((self.P.gamma-1.) * self.sol[0,1] / self.sol[0,0], (self.P.gamma-1.) * self.sol[-1,1] / self.sol[-1,0])
        str += "rhoe ee\t{0:10.3e}\t{1:10.3e}\n".format(self.P.peL/(self.P.gamma-1), self.P.peR/(self.P.gamma-1))
        # str += "    \t{0:10.3e}\t{1:10.3e}\n".format(self.sol[0,1], self.sol[-1,1])
        str += "-"*50 + "\n"
        str += "velocity of the traveling wave: sigma = {0:10.3e}\n".format(self.P.sigma)
        str += "-"*50 + "\n"
        str += "diffusion coefficient:       D = {:10.3e}\n".format(self.P.diffD)
        str += "thermal conductivity:   lambda = {:10.3e}\n".format(self.P.diffL)
        LD = self.P.diffD / abs(self.P.vhR - self.P.vhL)
        LL = (self.P.gamma-1.)/self.P.gamma*self.P.diffL/self.P.rhoeR / abs(self.P.vhR - self.P.vhL)
        str += "characteristic diffusion length: {:10.3e}\n".format(LD)
        str += "characteristic thermal length:   {:10.3e}\n".format(LL)
        str += "*"*80
        return str