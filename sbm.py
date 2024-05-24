import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def SBM(N, M, q0, q1):
    '''
    This function is designed to generate the Stochastic Block Model using networkx.
    input params (consistent with the project description):
    N (int): total number of nodes in the graph
    M (int): number of communities
    q0 (float): probability of intra-subgroup connection
    q1 (float): probability of inter-subgroup connection

    output:
    G: networkx graph object of the generated graph
    '''

    # Calculate the base size of each community and the number of extra nodes
    base_size = N // M
    extra_nodes = N % M

    # Initialize a list to hold the sizes of each community
    community_sizes = [base_size + 1 if i < extra_nodes else base_size for i in range(M)]

    # Create the list of sizes for the SBM generator
    sizes = community_sizes

    # Create the probability matrix for SBM
    p_matrix = np.full((M, M), q1)
    np.fill_diagonal(p_matrix, q0)

    # Generate the SBM using networkx
    G = nx.stochastic_block_model(sizes, p_matrix)

    return G

def visualize_graph(G):
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=50, edge_color='gray')
    plt.show()


if __name__ == "__main__":
    N = 100  # Total number of nodes
    M = 5    # Number of communities
    q0 = 0.8 # Intra-community connection probability
    q1 = 0.1 # Inter-community connection probability

    G = SBM(N, M, q0, q1)
    print(G)   # in latest networkx documentation, using print(G) works same as print(nx.info(G))
    visualize_graph(G)

