a1=1 
N=80 
K=40 
simul_instances_per_job=100000
na=20 
da=1 
L=60
alpha=.65
beta=.75
theta=.05
gamma_BO=.95
gamma_BS=.75
max_noise_amp=1
algorithm_option=2
try_max=20
e1=1
de=1
max_index=60
max_e=6

for err_bits in $(seq $e1 $de $max_e)
do
for index_var in $(seq $a1 $da $max_index) 
do 
qsub -N "clustered_neur_N_${N}${K}" -v N_var=80,K_var=40,L_var=60,max_inst=100000,err_bit=$err_bits,noise_amp=1,alpha0=0.65,beta0=0.75,theta0=0.05,gamma_BFO=0.85,gamma_BFS=0.82,algorithm_option=2,try_max=20,index=$index_var clustered_neural_recall.pbs 
done
done



