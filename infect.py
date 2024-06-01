import numpy as np
import random
import matplotlib.pyplot as plt
from sbm import SBM
import networkx as nx

def infect_step(G, p1, individuals, N):
    '''The function serves as the infection model for each day.
    input params (consistent with the project description):
    G (ndarray N*N): the adjacency matrix.
    p1: the probability each individual infects neighbours.
    individuals (ndarray): the current infection status of individuals.
    N (int): total number of individuals.

    output:
    individuals_updated (ndarray): the updated infection status of individuals.
    '''
    individuals_updated = individuals.copy()
    for i in range(N):
        if individuals[i] == 1:  # If the individual is infected
            neighbors = np.where(G[i] == 1)[0]  # Find all neighbors
            for neighbor in neighbors:
                if random.random() < p1:  # If neighbor is susceptible
                    individuals_updated[neighbor] = 1  # Infect the neighbor
    return individuals_updated

def infect(G, p0, p1, time_steps):
    '''The function serves as the infection model over multiple days.
    input params (consistent with the project description):
    G (ndarray N*N): the adjacency matrix.
    p0: the infection probability for initial status.
    p1: the probability each individual infects neighbours.
    time_steps: the number of days to simulate the infection spread.

    output:
    infection_progression (list of ndarray): the infection status of individuals over time.
    '''
    N = G.shape[0]
    individuals = np.random.choice([0, 1], size=N, p=[1-p0, p0])  # Initial infection status
    infection_progression = [individuals.copy()]  # Track the infection progression over time

    for _ in range(time_steps):
        individuals = infect_step(G, p1, individuals, N)
        infection_progression.append(individuals.copy())

    return infection_progression[time_steps]

def plot_infection_progression(infection_progression):
    days = len(infection_progression)
    infected_counts = [np.sum(day) for day in infection_progression]

    plt.figure(figsize=(10, 6))
    plt.plot(range(days), infected_counts, marker='o')
    plt.xlabel('Days')
    plt.ylabel('Number of Infected Individuals')
    plt.title('Infection Progression Over Time')
    plt.grid(True)
    plt.show()
