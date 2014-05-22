#=======================IMPORT THE NECESSARY LIBRARIES=========================
from brian import *
import time
import numpy as np
import sys
#==============================================================================


#================================INITIALIZATIONS===============================
n_exc = 320                             # The number of excitatory neurons in the output layer
n_inh = 80                              # The number of inhibitory neurons in the output layer
n = n_exc + n_inh                       # Total number of neurons in the output layer

synapse_delay = 0.0                     # The synaptic delay of ALL links in ms

connection_prob = 0.15                  # The probability of having a link from the input to output neurons in the second layer
p_minus = connection_prob * (float(n_inh)/float(n))
p_plus = connection_prob * (float(n_exc)/float(n))

input_stimulus_freq = 20000               # The frequency of spikes by the neurons in the input layer (in Hz)

frac_input_neurons = .3                # Fraction of neurons in the input layer that will be excited by a stimulus

no_samples_per_cascade = 3.0              # Number of samples that will be recorded
running_period = (no_samples_per_cascade/10.0)  # Total running time in seconds

no_cascades = 1000                    # Number of times we inject stimulus to the network
ensemble_size = 1                       # The number of random networks that will be generated

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

#======================READ THE NETWORK AND SPIKE TIMINGS======================
    
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
    
    
    #------------------------Read the Input Matrix-----------------------------
    file_name = "/Hesam/Academic/Network Tomography/Data/Graphs/We_cascades_%s.txt" %file_name_ending
    We = np.genfromtxt(file_name, dtype=None, delimiter='\t')
    
    file_name = "/Hesam/Academic/Network Tomography/Data/Graphs/Wi_cascades_%s.txt" %file_name_ending
    Wi = np.genfromtxt(file_name, dtype=None, delimiter='\t')
    
    W = np.vstack((We,Wi))
    #--------------------------------------------------------------------------
        

    #------------------------------Read and Sort Spikes------------------------
    file_name = "/Hesam/Academic/Network Tomography/Data/Spikes/S_times_cascades_%s.txt" %file_name_ending
    
    S_time_file = open(file_name,'r')
    S_times = np.fromfile(file_name, dtype=float, sep='\t')
    
    
    neuron_count = 0
    cascade_count = 0
    in_spikes = np.empty([n,no_cascades])
    out_spikes = np.empty([n,no_cascades])
    for l in range(0,len(S_times)-1):
        if (S_times[l] == -1):
            neuron_count = neuron_count + 1
        elif (S_times[l] == -2):
            neuron_count = 0
            cascade_count = cascade_count + 1
        else:
            if (S_times[l] < 0.00015):
                in_spikes[neuron_count][cascade_count] = 1
            else:
                out_spikes[neuron_count][cascade_count] = 1
    #--------------------------------------------------------------------------            
            
    #---------------------Check for Conflicts in Data--------------------------
    if (sum(multiply(in_spikes,out_spikes)) > 0):
        print('Error! Do something!')
        sys.exit()
    #--------------------------------------------------------------------------
        
#==============================================================================


#============================INFER THE CONNECTIONS=============================
    out_spikes_neg = -pow(-1,out_spikes)
    W_inferred = np.empty([n,n])
    for i in range(0,n-1):
        for j in range(0,n-1):
            W_inferred[i][j] = np.dot(in_spikes[i][:],out_spikes_neg[j][:])
    
#==============================================================================


#=============================TRANSFORM TO BINARY==============================
    W_binary = np.empty([n,n])
    for i in range(0,n-1):
        w_temp = W_inferred[:][i]
        ind = np.argsort(w_temp)
        w_temp = zeros(n)
        w_temp[ind[0:int(round(p_minus*n))]] = -1
        w_temp[ind[len(ind)-int(round(p_plus*n))+1:len(ind)]] = 1
        W_binary[:][i] = w_temp
        W_binary[i][i] = 0
#==============================================================================
    acc_plus = float(sum(multiply(W_binary>0,W>0)))/float(sum(W>0))
    acc_minus = float(sum(multiply(W_binary<0,W<0)))/float(sum(W<0))
    acc_zero = float(sum(multiply(W_binary==0,W==0)))/float(sum(W==0))
    file_name = "/Hesam/Academic/Network Tomography/Results/Accuracies/Acc_%s.txt" %file_name_ending
    acc_file = open(file_name,'a')
    acc_file.write("%f \t" % acc_plus)
    acc_file.write("%f \t" % acc_minus)
    acc_file.write("%f \t" % acc_zero)
    acc_file.write("\n")
    acc_file.close()

#=============================CALCULATE ACCURACY===============================


#==============================================================================
    