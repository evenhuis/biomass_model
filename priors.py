import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt


pr_range = {'O20'  : [100,400], 'DIC0' : [100,4000], 'TA0'  : [500,4000], \
            'N0'   : [600,1000], 'P0' : [10,50], \
            'mass0': [10,1000], 'pP0': [0.3,0.6],   'pC0': [0.0,0.3], \
            'P'    : [10,2000], 'R'    : [5,100],    'kla1' : [2,20]  ,    \
                                                   'kla2' : [5,200],     \
                   'fccm': [1.3,1.3], 'km1':[50,1000],'km2':[5,200], \
                   'dta':  [-100,100],\
                   'tauP': [1.5,15.], \
                   'tauR': [1.5,15.],\
                   'PQn' : [0.66,0.33], 'PQd' : [0.66,0.33],\
                   'NCd' : [0.,0.2 ],\
                   'PRd' : [-1,1], 'PRn' : [-1,1],\
                   'CRd' : [-1,1], 'CRn' : [-1,1],\
                   'LRd' : [-1,1], 'LRn' : [-1,1],\
            'O2_slope':[0.5,2.], \
            'O2_off':[-10,10], 'pH_off':[-0.5,0.5],  \
            'sg2_DA':[ 2,10],\
            'sg2_O2':[0.1,5], 'sg2_pH':[0.05,0.5]}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_gamma_pr( x, rng ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if( x <= 0. ):
        return -np.inf
    v0,v1=float(rng[0]),float(rng[1])
    mu =(v1+v0)/2.
    var=(v1-v0)**2/4.

    theta =   var/mu
    k     = mu**2/var
    return stats.gamma.logpdf( x,a=k, scale=theta )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_lognorm_pr( x, rng ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if( x <= 0 ):
        return -np.inf
    v0,v1 = float(rng[0]),float(rng[1])
    xb  = (v1+v0)/2.
    var = (v1-v0)**2/4.

    scale = xb/np.sqrt(1+var/xb**2)
    s2 = np.log(1+var/xb**2)
    return stats.lognorm.logpdf( x, np.sqrt(s2), scale=scale)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_norm_pr( x, rng ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    v0,v1 = float(rng[0]),float(rng[1])
    xb  = (v1+v0)/2.
    var = (v1-v0)**2/4.
    return stats.norm.logpdf( x-xb,  scale=var )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_smunif_pr( x, rng ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    loc = float(rng[0])
    scale = float(rng[1])-float(rng[0])
    alpha=0.90
    sg=(1.-alpha)/(alpha*np.sqrt(2*np.pi))*scale
    if x <= loc:
       return stats.norm.logpdf( x-loc,       scale=sg)+np.log(1-alpha)
    if loc+scale <= x :
       return stats.norm.logpdf( x-loc-scale, scale=sg)+np.log(1-alpha)
    return np.log(alpha/scale)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def prior_list( v, typ ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    global pr_range
    if(typ=='O20'     ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='DIC0'    ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='TA0'     ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='N0'      ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='P0'      ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='P'       ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='R'       ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='kla1'    ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='kla2'    ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='fccm'    ) : return stats.beta.logpdf( v, a=float(pr_range[typ][0]), b=float(pr_range[typ][1]))
    if(typ=='km1'     ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='km2'     ) : return log_gamma_pr( v, pr_range[typ] )
    if(typ=='tauP'    ) : return log_gamma_pr  ( v, pr_range[typ] )
    if(typ=='tauR'    ) : return log_gamma_pr  ( v, pr_range[typ] )
    if(typ=='fler'    ) : return log_lognorm_pr( v, pr_range[typ] )    
    if(typ=='PQd'     ) : return log_smunif_pr( v,  pr_range[typ]  )
    if(typ=='PQn'     ) : return log_smunif_pr( v,   pr_range[typ] )
    if(typ=='NCd'     ) : return log_norm_pr  ( v,   pr_range[typ] )
    if(typ=='O2_slope') : return log_gamma_pr( v, pr_range[typ ] )
    if(typ=='O2_off'  ) : return log_norm_pr ( v, pr_range[typ ] )
    if(typ=='pH_off'  ) : return log_norm_pr ( v, pr_range[typ ] )
    if(typ=='sg2_DA'  ) : return log_gamma_pr( v, pr_range[typ ] )
    if(typ=='sg2_O2'  ) : return log_gamma_pr( v, pr_range[typ ] )
    if(typ=='sg2_pH'  ) : return log_gamma_pr( v, pr_range[typ ] )
    return 1.0

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def plot_prior( typ, chain=None, nburn=150, rng=None ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if( rng is None ):
        xlo = 0.5*pr_range[typ][0]
        xhi = 1.5*pr_range[typ][1]
        nstep = 50
    else:
        xlo, xhi, nstep = rng

    if( chain is not None ):
        samps = np.reshape( chain[:,nburn:], [-1] )
        plt.hist( samps, bins=np.linspace( rng[0],rng[1],rng[2] ), normed=True )

    x  = np.linspace( xlo, xhi, nstep*10)
    px = list(map( lambda xx: prior_list(xx, typ ),x))
    plt.plot( x,np.exp(px))
    return

