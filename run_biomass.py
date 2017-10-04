import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

import sys
sys.path.append('../')
import os

from   biomass_model import py_interface_mod as pr
import priors
import scipy.optimize as op

import plot_biomass_model as pmr

from emcee_tools import helper_emcee as he

t0,t1 = 0,1

spline_control={}
spline_var    ={}

pr.py_setup_drivers("./driver.txt")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def interp_d( t, Ct, Cv, mode=0 ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''' interpolates or step inteps a spline'''
    if( mode==0 ):
        i = np.argmax( t<Ct)
        return Cv[i-1]
    if( mode==1 ):
        if( t<Ct[ 0]   ): return Cv[0]
        if(   Ct[-1]<t ): return Cv[-1]
        return np.interp( t, Ct, Cv)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def run_sim( theta, time  ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    t0=0
    ntime  = len(time)    ; t1=ntime +t0
    nlight = len(light_t) ; t2=nlight+t1
    nvalve = len(valve_t) ; t3=nvalve+t2
    nair   = len(air_t  ) ; t4=nair  +t3

    ta    = np.zeros(t4); mask =np.zeros(t4,dtype=bool)
    tc    = np.zeros(t4,int)
    ta[t0:t1] = time    ;tc[t0:t1]=0; mask[0:t4]=True
    ta[t1:t2] = light_t ;tc[t1:t2]=1;
    ta[t2:t3] = valve_t ;tc[t2:t3]=2; 
    ta[t3:t4] = air_t   ;tc[t3:t4]=3;

    isort = np.argsort( ta ) 
    ta    =   ta[isort] ; tc = tc[isort]
    mask  = mask[isort]

    tmask = np.logical_and( time[0]<=ta , ta<=time[-1] )
    #np.savetxt("timesteps.txt",np.column_stack((ta[tmask],tc[tmask])))

    yo = pr.sim_de( theta, ta[tmask]   )

    #return ta[np.logical_and(tmask,mask)],yo[:,mask[tmask]]
    return ta[tmask],yo

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def sim_on_obs( func, theta, x, x_o ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nsim    = len(x)                                # size of simulated variables
    [xu,iu] = np.unique( x_o, return_inverse=True ) # find unique obs only
    nobs_u  = len(xu)

    ta   = np.zeros( nobs_u+nsim )                  # big vector of time steps
    ta[0:nobs_u]   = xu
    ta[nobs_u: ]   = x

    mask = np.zeros( nobs_u+nsim, dtype=bool )    # T if an obs 
    mask[0:nobs_u] = True
   
    isort = np.argsort( ta )                      # sort time 
    ta   = ta[isort]
    mask = mask[isort]

    tt,ya = func( theta, ta )             # run the model
    return ya[:,mask][:,iu]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_likely( theta, x,O2_obs, pH_obs, DA_obs ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ta,y = run_sim( theta,x )
    yO2  = np.interp ( O2_obs[:,0],  ta, y[ 0] )
    yDIC = np.interp ( DA_obs[:,0],  ta, y[ 1] )
    yTA  = np.interp ( DA_obs[:,0],  ta, y[ 2] )
    ypH  = np.interp ( pH_obs[:,0],  ta, y[11] )

    ll =    np.sum(stats.norm.logpdf( yO2 -O2_obs[:,1], scale=theta[-2] )) \
           +np.sum(stats.norm.logpdf( ypH -pH_obs[:,1], scale=theta[-1] )) \
       +125*np.sum(stats.norm.logpdf( yDIC-DA_obs[:,1], scale=theta[-3]+0.50 )) \
       +125*np.sum(stats.norm.logpdf( yTA -DA_obs[:,2], scale=theta[-3]+0.50 )) 
    return ll

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def sg2_init( theta, x,O2_obs, pH_obs ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ta,y = run_sim( theta,x )
    yO2  = np.interp ( O2_obs[:,0],  ta, y[0] )
    ypH  = np.interp ( pH_obs[:,0],  ta, y[4] )
    return np.sqrt(np.average((yO2 - O2_obs[:,1])**2)), \
           np.sqrt(np.average((ypH - pH_obs[:,1])**2))


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_prior( theta ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    lp = 0
    for v,typ in zip(theta,theta_typ):
        #print(v,typ, prior_list( v,typ))
        lp += priors.prior_list( v,typ)
    return lp


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_prob( theta, *args ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    lprior = log_prior(theta)
    if( np.isfinite(lprior) ):
        llike = log_likely( theta, *args )
        if( np.isfinite(llike) ):
            return lprior + llike
        else:
            return -np.inf
    else:
        return -np.inf

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def nlogprob( theta, *args ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    nprob = -log_prob( theta, *args )
    if( np.isfinite(nprob) ):
        return  nprob
    else:
        return +1e+20


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def write_initial_entry(  theta, theta_typ, var ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    theta_typ = np.asarray(theta_typ)
    if( var in spline_control ):
        cont = spline_control[var]
        print( "{:10s}   {:d}  {:14.5g} {:14.5g}".format(var,len(cont),*list(map(float,priors.pr_range[var]) ) ))
        vals = theta[theta_typ==var]
        for c,v in zip(cont,vals):
            print("{:7.3f}  {:14.5g}".format(float(c),float(v)))
         
    else:
        print( "{:10s}   {:14.5g}  {:14.3g} {:14.3g}".format(var,*theta[theta_typ==var], *priors.pr_range[var] ) )
    return

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def write_initial_file(  theta, theta_typ ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    global t0, t1
    theta_typ = np.asarray(theta_typ)
    theta_typ_un, ind = np.unique(theta_typ,return_index=True)
    theta_typ_un = theta_typ_un[np.argsort(ind)]
    print("time_range {:6.3f} {:6.3f}".format(t0,t1))
    for var in theta_typ_un:
        write_initial_entry(  theta, theta_typ, var )   
    return

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def read_initial_entry( f, theta,th_typ ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    global spline_control, t0,t1
    variables  ='P R kla1 kla2 fccm km1 km2 dta tauP tauR fler PQd, PQn NCd PRd PRn CRd CRn LRd LRn'.split()
    var_singles='O20 DIC0 TA0 N0 P0 mass0 pPr0 pCb0 O2_slope O2_off pH_slope pH_off   sg2_DA sg2_O2 sg2_pH'.split()

    line = f.readline()
    if( line=="" ):
        return theta,th_typ,False
    line=line.split()

    var= line[0] 
    if( var=="time_range"):
        t0,t1 = float(line[1]),float(line[2])

    if( var in variables ):
        nvar =int(line[1])
        if(len(line)>2 ): 
            priors.pr_range[var]=line[2:]

        cont = np.zeros(nvar)
        vals = np.zeros(nvar)
        for i in range(nvar):
            th_typ.append(var)
            line=f.readline().split()
            cont[i],vals[i] = line[0:2]
        spline_control[var] = np.array(cont)
        spline_var    [var] = np.array(vals)

        pr.py_set_spline( var, cont )
        theta = np.append(theta,vals)

    if( var in var_singles ):
        if(len(line)>2 ):
            priors.pr_range[var]=list(map(float,line[2:]))
        theta = np.append(theta,float(line[1]))
        th_typ.append(var)
    return theta,th_typ,True

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def read_initial( fname ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    f = open(fname, 'r')
    theta   = np.zeros(0)
    theta_typ = []
    stat=True
    while stat:
        theta,theta_typ,stat = read_initial_entry(f,theta,theta_typ)
        if( not stat ):
            break
    return theta, theta_typ

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def get_data_trans( fname ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    f = open(fname, 'r')
    trans = np.zeros(4)
    trans[[0,2]]=1
    names='O2_slope O2_off pH_slope pH_off'.split()
    for line in f:
        parts = line.strip().split("=")
        key=parts[0].strip()
        try :
            i = names.index(key)
            trans[i] = float(parts[1])
        except:
            pass
    f.close()
    return trans



day=4

name=os.getcwd().split("/")[-1] +"/ day {} ".format(day)
theta,theta_typ = read_initial( "input_mod9_day{}.txt".format(day))

DA_obs = np.genfromtxt("DIC_Alk.txt")
O2_obs = np.genfromtxt("O2.txt")
pH_obs = np.genfromtxt("pH.txt")

trans = get_data_trans("data_trans.txt")
O2_obs[:,1] = trans[0]*(O2_obs[:,1]-210)+210 + trans[1]
pH_obs[:,1] = trans[2]*(pH_obs[:,1]-  7)+  7 + trans[3]


# read the valve_info

if( os.path.exists("valve_filt.txt") ):
    valve_t = np.genfromtxt("valve_filt.txt",usecols=[0])
else:
    valve_t = []

light_t = np.genfromtxt("I_filt.txt",usecols=[0])
air_t   = np.genfromtxt("air.txt",usecols=[0])

def mask_on_range( ts, rng, inv=False ):
    t0,t1 = rng
    mask = np.logical_and( t0< ts[:,0], ts[:,0]<t1 )
    if( inv ):
        mask = np.logical_not( mask)
    return ts[mask,:]

# time steps
time   = np.linspace(t0,t1,24*60)
O2_d = mask_on_range( O2_obs, [t0,t1] )
pH_d = mask_on_range( pH_obs, [t0,t1] )
DA_d = mask_on_range( DA_obs, [t0,t1] )

co2_treat = True 
if( not co2_treat ):
    # exclude the funny afternoon stuff
    tl0=day+16./24.
    tl1=day+27./24.
    O2_l = mask_on_range( O2_d, [tl0,tl1], inv=True )
    pH_l = mask_on_range( pH_d, [tl0,tl1], inv=True )
    DA_l = mask_on_range( DA_d, [tl0,tl1], inv=True )
else:
    O2_l =  O2_d
    pH_l =  pH_d
    DA_l =  DA_d



active = np.zeros( len(theta), dtype=bool )

t_mod, r_mod = run_sim( theta, time )
print("lp=",log_prob( theta, time,O2_l, pH_l, DA_l ))
pmr.plot_model_obs(  t_mod, r_mod, O2_d, pH_d, DA_d, name=name )
plt.savefig("plot_day{}.png".format(day))
plt.show()



# optimise just the varance terms
active[-3:] = True
ch = input('optimise? n')
if( not(ch in ["n","N"]) ):
	print("lp=",log_prob( theta, time,O2_l, pH_l, DA_l ))
	theta = he.optimise_ll_act( log_prob,  theta, active,time,O2_l, pH_l, DA_l,  nloop=10, nsub=5 )
	print("lp=",log_prob( theta, time,O2_l, pH_l, DA_l ))

# no activate them all
active[:] =True
#if( co2_treat ):
#    active[theta_typ.index("km1")]=False

simplex=None
while( not(ch in ["n","N"]) ):
    theta = he.optimise_ll_act( log_prob, theta, active, time,O2_l, pH_l, DA_l, live_plot=True, nloop=20, nsub=20)

   
    print(theta)
    max_LL = log_prob( theta, time,O2_l, pH_l, DA_l )
    print("lp=",max_LL )
    t_mod, r_mod = run_sim( theta, time )
    pmr.plot_model_obs(  t_mod, r_mod, O2_d, pH_d, DA_d, name=name )
    plt.savefig("plot_day{}.png".format(day))
    plt.show()

    ch = input('continue to optimise? n')



treat = "day{}_nm".format(day)
ch = input('sample chain? n')
if( not(ch in ["n","N"]) ):
    ch = input("Restart from chain? Y")
    if( ch in ["y","Y"] ):
        savefile="{}_chain.p".format(treat)
        if( os.path.isfile(savefile) ):
            probs,chain = he.load_chain(savefile)
            print("Chain loaded")
            p0 = chain[:,-1,:]
            samp = he.MCMC_restart( log_prob, p0, time, O2_l, pH_l, DA_l, nsteps=50, live_plot=True, max_LL=max_LL )
            fname = "{}_chain.p".format(treat)
            he.save_chain( fname, samp )
        else:
            print("chain file not found")
    else:
        samp = he.MCMC_all( log_prob, theta, time,  O2_l, pH_l, DA_l, nsteps=25, nwalker=200, live_plot=True, max_LL=max_LL  )
        fname = "{}_chain.p".format(treat)
        he.save_chain( fname, samp )

ch = input('Load the chain? Y')
if( ch in ["y","Y"] ):
    fname = "{}_chain.p".format(treat)
    [prob,chain] = he.load_chain( fname )
