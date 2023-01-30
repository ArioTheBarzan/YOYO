import random
import numpy as np
import statistics
import matplotlib.pyplot as plt
import time

start_time = time.time()

#def shit_test(shit_list):
#    for shit in range(len(shit_list)):
#        for crap in range(len(shit_list)):
 #           if crap != shit and shit_list[shit] == shit_list[crap]:
#                print("errooooooor")
#                print(shit_list)
 #               exit()

def En_of_S(jay, state, length, bond_vert, bond_hor):
    energy = 0
    delE = -0.5 * jay
    number = length * length
    for i in range(number):
        energy += delE * bond_vert[i] * state[i] * state[(i - length) % number]
        energy += delE * bond_hor[i] * state[i] * state[(i // length) * length + (i - 1) % length]

    return energy


L = 32
N = L * L

J = 1.0  # strength of mag field
T = 2.2  # temperature
p = 1.0 - np.exp(-2.0 / T)  # probability of wolff
p2 = 0.5  # prob of bond frustration

plotp = []
plot_corr = []

center = 0
rx = 16
ry = 16
corr_center_to_r = 0

nbr = {i: ((i // L) * L + (i + 1) % L, (i + L) % N, (i // L) * L + (i - 1) % L, (i - L) % N) for i in
       range(N)}  # defines a neighbourhood\

l_row = 0
l_coloumn = 0

tau_row = [1 for k1 in range(N)]
tau_coloumn = [1 for k2 in range(N)]

nsteps = 10000
ntsteps = 200
tot_pstep = 20
tot_nrun = 20
p_inc = p2 / tot_pstep



for pstep in range(tot_pstep):
    print("------------------------------------------------------------------------")
    print((time.time() - start_time) // 3600, "and", ((time.time() - start_time) % 3600)/60)
    p2 -= p_inc
    tau_average_corr = 0
    print("p2 is ", p2, " while step number is ", pstep)
    for nrun in range(tot_nrun):
        S = [random.choice([1, -1]) for k in range(N)] # defines initial conditions for the lattice
        tau_row = [1 for k1 in range(N)]
        tau_coloumn = [1 for k2 in range(N)]

        for x in range(N):  # defines tau structure for each step
            if random.uniform(0.0, 1.0) < p2:  # with probability p2 assigns frustrated bonds
                tau_coloumn[x] = -1
            if random.uniform(0.0, 1.0) < p2:
                tau_row[x] = -1
        therm_average_corr = 0  # thermal average of the spin-spin correlation for a specific tau structure
        Z = 0
        for step in range(nsteps):
            E = En_of_S(1.0, S, L, tau_coloumn, tau_row)
            therm_average_corr += S[0] * S[(6 + 6 * L) % N] * np.exp(-1 * E / T)
            Z += np.exp((-1 * E) / T)
            k = random.randint(0, N - 1)
            Pocket, Cluster = [k], [k]
            while Pocket != []:
                j = random.choice(Pocket)
                for l in nbr[j]:
                    if l not in Cluster:

                        if random.uniform(0.0, 1.0) < p:
                            if l == (j - L) % N:
                                l_coloumn = j

                                if tau_coloumn[l_coloumn] == -1 and S[l] == -S[j] and l not in Cluster:
                                    Pocket.append(l)
                                    Cluster.append(l)

                                if tau_coloumn[l_coloumn] == 1 and S[l] == S[j] and l not in Cluster:
                                    Pocket.append(l)
                                    Cluster.append(l)

                            if l == (j // L) * L + (j - 1) % L:
                                l_row = j

                                if tau_row[l_row] == 1 and S[l] == S[j] and l not in Cluster:
                                    Pocket.append(l)
                                    Cluster.append(l)

                                if tau_row[l_row] == -1 and S[l] == -S[j] and l not in Cluster:
                                    Pocket.append(l)
                                    Cluster.append(l)

                            if l == (j // L) * L + (j + 1) % L:
                                l_row = l

                                if tau_row[l_row] == 1 and S[l] == S[j] and l not in Cluster:
                                    Pocket.append(l)
                                    Cluster.append(l)

                                if tau_row[l_row] == -1 and S[l] == -S[j] and l not in Cluster:
                                    Pocket.append(l)
                                    Cluster.append(l)

                            if l == (j + L) % N:
                                l_coloumn = l

                                if tau_coloumn[l_coloumn] == -1 and S[l] == -S[j] and l not in Cluster:
                                    Pocket.append(l)
                                    Cluster.append(l)

                                if tau_coloumn[l_coloumn] == 1 and S[l] == S[j] and l not in Cluster:
                                    Pocket.append(l)
                                    Cluster.append(l)

                Pocket.remove(j)
            for j in Cluster:
                S[j] *= -1
        therm_average_corr /= Z
        tau_average_corr += therm_average_corr / tot_nrun  # tau structure average of the thermal average of the spin-spin correlation for a set p2



    plotp.append(p2)
    plot_corr.append(tau_average_corr)

fig, ax = plt.subplots()
plt.plot(plotp, plot_corr, color='r', label='correlation')
plt.title("T = 2.2")
plt.xlabel("Probability")
plt.ylabel("Correlation to 16 acroos and 16 down")
plt.legend()
plt.show()
plt.savefig("output-1.jpg")


print("thermal average is ", therm_average_corr)

end_time = time.time()

print("total runtime was : ", end_time - start_time)