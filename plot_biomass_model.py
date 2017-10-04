import matplotlib.pyplot as plt
import numpy as np

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def plot_model_obs(   time, mod_res, O2_o=None, pH_o=None, DA_o=None, sim_chain=None, name=None ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    time_col  ="#189ded"
    o2_col  ="#cf386f"
    dic_col ="#757575"
    co2_col ='#f75d0b'
    hco3_col="#449d26"
    co3_col ="#427548"

    # run the model
    fig,axs= plt.subplots(4, 2, sharex=True,figsize=(5, 9))

    if( sim_chain is not None ):
        [ns,nv,nt] = np.shape(sim_chain)
        ntime = np.size(time)
        print(nt,ntime)
        for i in range(ns):
            axs[0,0].plot( time, sim_chain[i, 0,:], color='red',alpha=0.2)
            axs[1,0].plot( time, sim_chain[i,11,:], color='red',alpha=0.2)

            axs[2,0].plot( time, sim_chain[i, 1,:], color=dic_col,alpha=0.2)
            axs[2,0].plot( time, sim_chain[i, 2,:], color=time_col,alpha=0.2)
            axs[2,0].plot( time, sim_chain[i, 5,:], color=co2_col,alpha=0.2)
            axs[2,0].plot( time, sim_chain[i, 6,:], color=hco3_col,alpha=0.2)

            axs[3,0].plot( time, sim_chain[i, 9,:]/24., color="blue",alpha=0.2)
            axs[3,0].plot( time, sim_chain[i,10,:]/24., color="red",alpha=0.2)

    axs[0,1].set_ylabel("DW (mg/L)")
    axs[0,1].set_ylim(0,250)
    axs[0,1].plot( time, mod_res[18]/1e3)
    axs[0,1].plot( time, mod_res[19]/1e3)
    axs[0,1].plot( time, mod_res[20]/1e3)
    axs[0,1].plot( time, (mod_res[18]+mod_res[20]+mod_res[19])/1e3)
    axs[0,1].grid()

    axs[2,1].set_ylabel(r"Conc ($\mu$M)")
    axs[2,1].set_ylim(0,1000)
    axs[2,1].plot( time, mod_res[3],color='red')
    axt = axs[2,1].twinx()
    axt.plot( time, mod_res[4],color='green')
    axt.set_ylim(0,50)
    axs[2,1].grid()    



    axs[0,0].set_xlim(np.min(time),np.max(time))
    axs[0,0].axhline( 210 , ls='--',color='red')
    axs[0,0].set_ylabel("O$_2$ / $\mu$M")
    if( O2_o is not None ) : 
        axs[0,0].plot( O2_o[:,0],O2_o[:,1] )
    axs[0,0].plot( time, mod_res[0,:], '-', color='red' )
    axs[0,0].grid()
    #axs[0].plot( time, mod_res[0,:],'o',color='red', ms=1.0, ls=None)

    axs[1,0].set_ylabel("pH")
    axs[1,0].set_ylim(7.5,10.5)
    axs[1,0].grid()
    if( pH_o is not None ): 
        axs[1,0].plot( pH_o[:,0],pH_o[:,1] )
              
    axs[1,0].plot( time, mod_res[11,:], color='red' )

    axs[2,0].set_ylabel("conc / $\mu$M")
    if( DA_o is not None ):
        axs[2,0].plot( DA_o[:,0], DA_o[:,1],  'o',color=dic_col)
        axs[2,0].plot( DA_o[:,0], DA_o[:,2],  'o',color=time_col)
        axs[1,0].plot( DA_o[:,0], DA_o[:,5],  'o',color='blue')
        axs[0,0].plot( DA_o[:,0], DA_o[:,3],  'o',color='blue') 
    axs[2,0].plot( time, mod_res[1,:], color=dic_col ,label="DIC" )
    axs[2,0].plot( time, mod_res[2,:], color=time_col  ,label="TA" )
    axs[2,0].plot( time, mod_res[5,:], color=co2_col ,label="CO$_2$" )
    axs[2,0].plot( time, mod_res[6,:], color=hco3_col,label="HCO$_3$" )
    axs[2,0].plot( time, mod_res[7,:], color=co3_col ,label="CO$_3$" )
    axs[2,0].grid()
    #axs[2,0].legend()

    axs[3,0].set_ylabel("rate / $\mu$M/hr")
    axs[3,0].set_ylim(-100,500)
    axs[3,0].grid()
    axs[3,0].plot( time, mod_res[-2,:]/24., label="P" )
    axs[3,0].plot( time, mod_res[-1,:]/24., label="R" )

    plt.tight_layout()
    if( name is not None) :
        plt.suptitle(name)
    return
