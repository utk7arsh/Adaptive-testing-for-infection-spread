from infect import infect
import numpy as np
import random
import networkx as nx
from sbm import SBM
import matplotlib.pyplot as plt



def naive_testing(s):
    num = 0
    stages = 1
    new_s = np.zeros(len(s))
    for i in range(len(s)):
        if s[i] == 1:
            new_s[i] = 1
        num += 1
    return num, stages, new_s


def optimized_naive_testing(s):
    # Fix test count
    num = 0
    stages = 1
    positive_tests = 0
    negative_tests = 0

    proportion_of_positive_tests = sum(s)/len(s)
    num+=1
    
    if(proportion_of_positive_tests <= 0.5):
      new_s = np.zeros(len(s))
      for i in range(len(s)):
          if s[i] == 1:
            new_s[i] = 1
            positive_tests += 1
          num+=1
          if positive_tests==sum(s): break
    else:
      new_s = np.ones(len(s))
      for i in range(len(s)):
          if s[i] == 0:
            new_s[i] = 0
            negative_tests += 1
          num+=1
          if negative_tests==(len(s)-sum(s)): break

    return num, stages, new_s

# binary spliting
def binary_splitting_round(s):
    # s: np.array the infectious status & test status
    num = 0
    flag = sum(s[:,0])>0
    assert flag
    stages = 0
    if len(s[:,0])==1:
        s[0,1] = s[0,0]
        return num,s,stages
    
    B1, B2 = np.array_split(s.copy(), 2,axis=0)
    flag = sum(B1[:,0])>0
    num+=1
    stages += 1
    
    if flag:
        n,stmp,stage = binary_splitting_round(B1)
        s[:len(B1),1] = stmp[:,1]
    else:
        s[:len(B1),1] = 0
        n,stmp,stage = binary_splitting_round(B2)
        s[len(B1):,1] = stmp[:,1]
    num += n
    stages += stage
    return num,s,stages 

def binary_splitting(s):
    # modified bs
    # s: 1-d array the infectious status
    st = np.zeros((len(s),2))
    st[:,0] = s
    st[:,1] = np.nan
    nums = 0
    count = sum(np.isnan(st[:,1]))
    stages = 0
    # the undetermined people
    while count!=0:
        mask = np.isnan(st[:,1])
        flag = sum(st[mask,0]>0)>0
        nums += 1
        stages+=1
        if not flag:
            st[mask,1] = 0
        else:
            n,stmp,stage = binary_splitting_round(st[mask,:])
            st[mask,1] = stmp[:,1]
            nums += n
            stages += stage
        count = sum(np.isnan(st[:,1]))
        
    assert sum(st[:,0]!=st[:,1])==0
    return nums,stages, st[:,1]

# diag
def diagalg_iter(s):
    # s(np.array): binary string of infection status
    k = int(np.log2(len(s)))
    l = int(2**(k-1))
    lp = 0
    p = np.zeros(k+1)
    group = dict()
    num = np.ones(k+1,dtype=np.int32)
    for i in range(k):
        p[i] = sum(s[lp:lp+l])>0
        group[i] = s[lp:lp+l]
        num[i] = l
        lp+=l
        l = l//2

    p[-1] = s[-1]
    group[k] = np.array([s[-1]])
    # p(array): pattern
    # group(dict): indicate the group information
    # num(array): the group size
    return p.astype(np.int32), group,num


def diag_splitting(s):
    # s(np.array): binary string of infection status
    num_tests = 0
    stages = 0
    pattern, group, nums = diagalg_iter(s)
    stages +=1
    num_tests += len(pattern)
    indices = np.where(pattern == 1)[0]
    flag = 0
    for i in indices:
        if nums[i]>1:
            num_test,stage = diag_splitting(group[i])
            num_tests += num_test
            if not flag:
                stages+=stage
                flag = 1
    return num_tests,stages

def Qtesting1_iter(s):
    # s(np.array): binary string of infection status
    # Ex: [1,0,0,1,0,0,1,0]
    # k = 3 (iteration number)
    k = int(np.log2(len(s)))
    # l = 4
    l = int(2**(k-1))
    lp = 0
    # p = [0,0,0,0] (num_of_groups)
    p = np.zeros(k+1)
    group = dict()
    # num = [1,1,1,1]
    num = np.ones(k+1,dtype=np.int32)
    for i in range(k):
        p[i] = sum(s[lp:lp+l])>0
        group[i] = s[lp:lp+l]
        num[i] = l
        lp+=l
        l = l//2

    p[-1] = s[-1]
    group[k] = np.array([s[-1]])
    # p(array): pattern
    # group(dict): indicate the group information
    # num(array): the group size
    return p.astype(np.int32), group,num

def Qtesting1(s):
    # This is technically a test but I account for it within optimized_naive_testing (instead of passing a parameter)
    if sum(s)/len(s) > 0.2:
        num, stages, _ = optimized_naive_testing(s)
        return num, stages

    # s(np.array): binary string of infection status
    num_tests = 0
    stages = 0
    pattern, group, nums = Qtesting1_iter(s)
    stages +=1
    # pattern is the test_result for each group
    num_tests += len(pattern)

    # All groups which tested positive
    indices = np.where(pattern == 1)[0]
    flag = 0
    for i in indices:
        if nums[i]>1:
            num_test,stage = Qtesting1(group[i])
            num_tests += num_test
            if not flag:
                stages+=stage
                flag = 1
    return num_tests,stages


def Qtesting2(s):
    '''
    s(np.array): binary string of infection status
    '''
    num_tests = 0
    stages = 0
    ###################################################
    '''your code here'''

    ###################################################



    return num_tests,stages



def Qtesting1_comm_aware(s, communities):
    '''
    s(np.array): binary string of infection status
    communities(list): the community information
    '''
    num_tests = 0
    stages = 0

    return num_tests, stages

# Main function to compare the performance with binary, diagonal, and naive testing
def main():
    infection_probabilities = np.linspace(0.01, 0.3, 30)
    num_tests_naive = []
    num_tests_binary = []
    num_tests_diag = []
    num_tests_q1_comm_aware = []
    stages_naive = []
    stages_binary = []
    stages_diag = []
    stages_q1_comm_aware = []

    N = 50  # Total number of nodes
    M = 5   # Number of communities
    q0 = 0.8  # Intra-community connection probability
    q1 = 0.1  # Inter-community connection probability
    G = SBM(N, M, q0, q1)
    adj_matrix = nx.to_numpy_array(G)

    communities = []
    for community in nx.connected_components(G):
        communities.append(list(community))

    for p in infection_probabilities:
        array = np.random.choice([0, 1], size=N, p=[1 - p, p])

        num, stages, _ = naive_testing(array)
        num_tests_naive.append(num)
        stages_naive.append(stages)

        num, stages, _ = binary_splitting(array)
        num_tests_binary.append(num)
        stages_binary.append(stages)

        num, stages = diag_splitting(array)
        num_tests_diag.append(num)
        stages_diag.append(stages)

        num, stages = Qtesting1_comm_aware(array, communities)
        num_tests_q1_comm_aware.append(num)
        stages_q1_comm_aware.append(stages)

    plt.figure(figsize=(16, 8))

    # Plot Number of Tests
    plt.subplot(1, 2, 1)
    plt.plot(infection_probabilities, num_tests_naive, label="Naive Testing", color='blue')
    plt.plot(infection_probabilities, num_tests_binary, label="Binary Splitting", color='orange')
    plt.plot(infection_probabilities, num_tests_diag, label="Diagonal Splitting", color='green')
    plt.plot(infection_probabilities, num_tests_q1_comm_aware, label="Q1 Testing (Comm Aware)", color='red')
    plt.xlabel("Infection Probability")
    plt.ylabel("Number of Tests")
    plt.title("Number of Tests vs Infection Probability")
    plt.legend()
    plt.grid(True)

    # Plot Number of Stages
    plt.subplot(1, 2, 2)
    plt.plot(infection_probabilities, stages_naive, label="Naive Testing", color='blue')
    plt.plot(infection_probabilities, stages_binary, label="Binary Splitting", color='orange')
    plt.plot(infection_probabilities, stages_diag, label="Diagonal Splitting", color='green')
    plt.plot(infection_probabilities, stages_q1_comm_aware, label="Q1 Testing (Comm Aware)", color='red')
    plt.xlabel("Infection Probability")
    plt.ylabel("Number of Stages")
    plt.title("Number of Stages vs Infection Probability")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()


def Qtesting2_comm_aware(s,communities):
    '''
    s(np.array): binary string of infection status
    communities(list): the community information
    '''
    num_tests = 0
    stages = 0
    ###################################################
    '''your code here'''

    ###################################################



    return num_tests,stages