#=======================IMPORT THE NECESSARY LIBRARIES=========================
from brian import *
import time
import numpy as np 
#==============================================================================


#================================INITIALIZATIONS===============================
n_exc = 320                             # The number of excitatory neurons in the output layer
n_inh = 80                              # The number of inhibitory neurons in the output layer
n = n_exc + n_inh                       # Total number of neurons in the output layer

synapse_delay = 0.0                       # The synaptic delay of ALL links in ms

connection_prob = 0.15                  # The probability of having a link from the input to output neurons in the second layer


input_stimulus_freq = 20000               # The frequency of spikes by the neurons in the input layer (in Hz)

frac_input_neurons = .3                # Fraction of neurons in the input layer that will be excited by a stimulus

no_samples_per_cascade = 3.0              # Number of samples that will be recorded
running_period = (no_samples_per_cascade/10.0)  # Total running time in seconds

no_cascades = 400                    # Number of times we inject stimulus to the network
ensemble_size = 5                       # The number of random networks that will be generated

#----------------------------------Neural Model--------------------------------
tau=10*ms
tau_e=2*ms # AMPA synapse
eqs='''
dv/dt=(I-v)/tau : volt
dI/dt=-I/tau_e : volt
'''
#------------------------------------------------------------------------------


def myrates(t):
    rates=zeros(n)*Hz    
    if t < 0.1 * ms:
        input_index = floor(n*rand(round(n*frac_input_neurons)))
        input_index = input_index.astype(int)
        rates[input_index]=ones(round(n*frac_input_neurons))*input_stimulus_freq *Hz
    return rates
    
#==============================================================================




for ensemble_count in range(0,ensemble_size):

#============================GENERATE THE NETWORK==============================
    
    
    
    #----------------------Initialize the Main Layer---------------------------
    main_network=NeuronGroup(n,model=eqs,threshold=10*mV,reset=0*mV)

    Pe = main_network.subgroup(n_exc)
    Pi = main_network.subgroup(n_inh)
    Ce = Connection(Pe, main_network, weight=1*mV, sparseness=connection_prob,delay = synapse_delay*ms)
    Ci = Connection(Pi, main_network, weight=-1*mV, sparseness=connection_prob,delay = synapse_delay*ms)
        
        
    #--------------------------------------------------------------------------
    
    
    #--------------Transform Connections to Weighted Matrices------------------
    #Wf = input_connections.W.todense()
    We = Ce.W.todense()
    Wi = Ci.W.todense()
    #--------------------------------------------------------------------------
    
    
    #----------------------Construct Prpoper File Names------------------------
    file_name_ending = "n_exc_%s" %str(int(n_exc))
    file_name_ending = file_name_ending + "_n_inh_%s" %str(int(n_inh))    
    file_name_ending = file_name_ending + "_p_%s" %str(connection_prob)
    file_name_ending = file_name_ending + "_r_%s" %str(frac_input_neurons)
    file_name_ending = file_name_ending + "_f_%s" %str(input_stimulus_freq)
    file_name_ending = file_name_ending + "_d_%s" %str(synapse_delay)
    file_name_ending = file_name_ending + "_T_%s" %str(no_cascades)    
    file_name_ending = file_name_ending + "_%s" %str(ensemble_count)
    #--------------------------------------------------------------------------
    
    #----------------Save Connectivity Matrices to the File--------------------
    #np.savetxt("/Hesam/Academic/Network Tomography/Data/Graphs/Wf_cascades_%s.txt" %file_name_ending ,Wf,'%1.4f',delimiter='\t',newline='\n')
    np.savetxt("/Hesam/Academic/Network Tomography/Data/Graphs/We_cascades_%s.txt" %file_name_ending,We,'%1.4f',delimiter='\t',newline='\n')
    np.savetxt("/Hesam/Academic/Network Tomography/Data/Graphs/Wi_cascades_%s.txt" %file_name_ending,Wi,'%1.4f',delimiter='\t',newline='\n')
    #--------------------------------------------------------------------------
 
    #-------------------Run the Network and Record Spikes----------------------
    S_time_file = open("/Hesam/Academic/Network Tomography/Data/Spikes/S_times_cascades_%s.txt" %file_name_ending,'w')
    for cascade_count in range(0,no_cascades):
        
        #-------------------Initialize the Input Stimulus----------------------
        
        inputer_dummy_layer=PoissonGroup(n,myrates)
        input_connections=Connection(inputer_dummy_layer,main_network,weight=lambda i,j:(1-abs(sign(i-j))),delay = 0*ms)
        #----------------------------------------------------------------------
        
        
        M_l1 = SpikeMonitor(inputer_dummy_layer)
        M_l2 = SpikeMonitor(main_network)
        M_l1.source.clock.reinit()                  # Reset the spikemonitor's clock so that for the next random network, everything starts from t=0
        M_l2.source.clock.reinit()                  # Reset the spikemonitor's clock so that for the next random network, everything starts from t=0
        run(running_period * ms)        
    
        print M_l1.nspikes, "spikes in input layer"
        print M_l2.nspikes, "spikes in output layer"
        #----------------------------------------------------------------------
    
        #--------------------Save Spike Times to the File----------------------        
        
        
        SS = M_l2.spiketimes
        for l in range(0,len(SS)-1):
            item = SS[l]
            
            for j in range(0,len(item)):
                if (item[j]<0.00011):
                    S_time_file.write("%f \t" % item[j])
                if (item[j]>0.00011):
                    S_time_file.write("%f \t" % item[j])
                    break
            S_time_file.write("-1 \n")
    
        if (cascade_count<no_cascades-1):
            S_time_file.write("-2 \n")
        #----------------------------------------------------------------------

            M_l1.reinit()                               # Reset the spikemonitor so that for the next random network, everything starts from t=0
            M_l2.reinit()                               # Reset the spikemonitor so that for the next random network, everything starts from t=0
            M_l1.source.reinit()
            M_l2.source.reinit()
            M_l1.source.reset()            
            M_l2.source.reset()
    S_time_file.close()