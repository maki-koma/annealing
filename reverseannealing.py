from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import numpy as np
token = ''
endpoint = 'https://cloud.dwavesys.com/sapi/'

A=50
B=30
C=10
D=10
#T,時刻, K, qubit数, N,台数
T=int(input("入力"))
K=2*T
N=int(input("入力"))

x=np.array([1, 1, 0, 0, 1, 1, 0, 0, 1, 0])

import matplotlib.pyplot as plt
t = (0,5,15,20)
s = (1.0,0.2,0.2,1.0)
plt.plot(t,s)
plt.show()

schedule = list(zip(t,s))
print(schedule)

while D <160:
    H = np.zeros(K*K).reshape(K,K)
    for i in range(T-1):
        H[(i,i+1)] = A
    for i in range(T, 2*T-1):
        H[(i,i+1)] = A
   
    for i  in range(T-2):
        H[(i,i+2)] = B
    for i  in range(T, 2*T-2):
        H[(i,i+2)] = B
   
    for i  in range(T):
        H[(i,i+5)] = C
       
    for i in range(K):
        for j in range(K):
            if i < j:
                H[(i,j)] = H[(i,j)]+2*D
            if i == j:
                H[(i,j)] = H[(i,j)]-2*N*D+D
                
    dw_sampler = DWaveSampler(solver='Advantage_system1.1', token=token)
    sampler = EmbeddingComposite(dw_sampler)
    sampleset = sampler.sample_qubo(H,num_reads=100, anneal_schedule=schedule, initial_state = x)
    for sample, energy, num_occurrences, chain_break_fraction in list(sampleset.data()):
        print("Dの値は", D, sample, "Energy: ", energy, "Occurrences: ", num_occurrences)
    D += 10
