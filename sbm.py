import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def SBM(N, M, q0, q1):
    """
    This function is designed to generate the Stochastic Block Model.

    Parameters:
    N (int): The number of individuals
    M (int): The number of subgroups
    q0 (float): Probability of intra-subgroup connection
    q1 (float): Probability of inter-subgroup connection

    Returns:
    np.ndarray: Adjacency matrix of the generated graph
    """

    G = np.zeros((N, N))

    # Generate node indices and split them into subgroups...
    nodes = np.arange(N)
    communities = np.array_split(nodes, M)

    # Iterate over communities to create connections.
    for i in range(M):
        comm1 = communities[i]
        for j in range(i, M):
            comm2 = communities[j]

            if i == j:
                # Generate intra-community connections
                if q0 > 0:
                    intra_block = np.random.choice([0, 1], size=(comm1.size, comm1.size), p=[1 - q0, q0])
                    G[np.ix_(comm1, comm1)] = np.triu(intra_block, 1) + np.triu(intra_block, 1).T
            else:
                # Generate inter-community connections
                if q1 > 0:
                    inter_block = np.random.choice([0, 1], size=(comm1.size, comm2.size), p=[1 - q1, q1])
                    G[np.ix_(comm1, comm2)] = inter_block
                    G[np.ix_(comm2, comm1)] = inter_block.T

    return G


def visualize_graph(G):
    """
    Function to visualize the adjacency matrix of a graph generated by the Stochastic Block Model.

    Parameters:
    G (np.ndarray): Adjacency matrix of the graph.

    Returns:
    None
    """
    graph = nx.from_numpy_array(G)

    # Draw the graph
    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(graph, seed=42)  # Position nodes using a spring layout
    nx.draw(graph, pos, with_labels=True, node_color='skyblue', node_size=500, edge_color='gray', linewidths=1, font_size=15)
    plt.title("Visualization of Stochastic Block Model (SBM) Graph")
    plt.show()


if __name__ == "__main__":
    N = 50  # Total number of nodes
    M = 5    # Number of communities
    q0 = 1 # Intra-community connection probability
    q1 = 0 # Inter-community connection probability

    G = SBM(N, M, q0, q1)
    print(G)   # in latest networkx documentation, using print(G) works same as print(nx.info(G))
    visualize_graph(G)