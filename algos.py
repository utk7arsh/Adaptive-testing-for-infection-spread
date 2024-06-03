from infect import infect
import numpy as np
import random
import networkx as nx
from sbm import SBM
import matplotlib.pyplot as plt


def individual_testing(s):
    num = len(s)
    stages = 1
    return num, stages


def optimized_individual_testing(s):
    """
    Aim: This is an optimization over individual testing by leveraging T1 test information
    Features: Search Pruning, 0-1 Symmetry
    Performance: Approx 30% improvement over individual testing when infection rate is extremely 
                low or high (<0.1 or 0.9>)
    """
    num = 0
    stages = 1
    total_infected = sum(s)
    num+=1
    proportion_of_positive_tests = total_infected/len(s)

    # To avoid hassle of individual testing
    if proportion_of_positive_tests==0 or proportion_of_positive_tests==1: 
        return num, stages

    stages += 1
    positive_tests = 0
    negative_tests = 0

    # Based on proportion, we either look for ones or zeros
    if(proportion_of_positive_tests <= 0.5):
      for i in range(len(s)):
          if s[i] == 1:
            positive_tests += 1
          num+=1
          if positive_tests==sum(s): break

    else:
      for i in range(len(s)):
          if s[i] == 0:
            negative_tests += 1
          num+=1
          if negative_tests==(len(s)-sum(s)): break

    return num, stages



# binary spliting
def binary_splitting_round(s):
    # s: np.array the infectious status & test status
    # Number of tests performed
    num = 0
    # num infected > 1?
    flag = sum(s[:,0])>0
    assert flag
    # num stages till now
    stages = 0
    # if there's one person, set their test status to infection status
    if len(s[:,0])==1:
        s[0,1] = s[0,0]
        return num,s,stages

    # Split array in half
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



def parallel_binary_splitting(s):

    """
    Aim: Reduce repetitive overhead of binary splitting and leverage T1 test information
    Features: Parallel Searching, Search Pruning, 0-1 Symmetry
    Performance: Way beter than binary splitting, on average worse than diagonal splitting
    [CONSIDER ADDING THRESHOLDING BASE CASE - idea discarded after poor results]
    """

    num, stages = 0,1
    if len(s)==0: return num, stages

    num+=1
    num_infected = sum(s)
    # Base Case: Prune len(s)
    if num_infected == 0 or num_infected==len(s): return num, stages
    # if num_infected > 0.8: return optimized_individual_testing(s)

    mid = len(s)//2
    left_half = s[:mid].copy()
    right_half = s[mid:].copy()
    
    num_left, stages_left = parallel_binary_splitting(left_half)
    num_right, stages_right = parallel_binary_splitting(right_half)
    
    num += num_left + num_right
    stages += max(stages_left,stages_right)

    return num, stages



# diag
def diagalg_iter(s):
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


def group_division_testing(s, threshold = 0.15):
    """
    Aim: Divide the population into groups of size s before testing. 
    Features: Group size determined by infection rate, 
            Switch to parallel_binary_splitting if infection rate is too low
            Switch to individual testing (optimized) if infection rate is too high
    Performance: Significantly better than diagonal testing for all initial probabilities!
    """

    num, stages = 0, 0
    if len(s) == 0: 
        return num, stages

    num_infected = sum(s)
    num+=1
    stages += 1

    if num_infected == 0 or num_infected==len(s): return num, stages

    if num_infected/len(s) < 0.001 or num_infected/len(s) > 0.999:
        bin_num, bin_stages = parallel_binary_splitting(s)
        num += bin_num
        stages += bin_stages
        return num, stages
    
    # Lesser the proportion of infected, larger the group size?
    if num_infected/len(s) < threshold or num_infected/len(s) > (1-threshold):
        group_size = max(1, len(s) // 32)
    else:
        group_size = max(1, len(s) // 64)
    groups = [s[i:i + group_size] for i in range(0, len(s), group_size)]
    
    # Apply conditional parallel binary splitting on each group
    for group in groups:
        # if len(group) > 1:
            num_sub, stages_sub = Qtesting1_bin(group)
            num += num_sub
            stages = max(stages, stages_sub)
    
    return num, stages

def perform_T2_test(s):
    num_infected = sum(s)
    # Key:
    # 0 --> [0,0]
    # 1 --> [1,2)
    # 3 --> [2,4)
    # 6 --> [4,8)
    # _ --> [8, len(s))
    if num_infected == 0: return 0
    elif num_infected < 2: return 1
    elif num_infected < 4: return 3
    elif num_infected < 8: return 6
    else: return (8+len(s))/2


def parallel_binary_splitting_T2(s):
    """
    Aim: Adapt parallel binary splitting to use only T2 tests
    """

    num, stages = 0,1
    if len(s)==0: return num, stages

    num+=1
    num_infected = perform_T2_test(s)
    if num_infected == 0 or len(s)==1: return num, stages

    mid = len(s)//2
    left_half = s[:mid].copy()
    right_half = s[mid:].copy()
    
    num_left, stages_left = parallel_binary_splitting_T2(left_half)
    num_right, stages_right = parallel_binary_splitting_T2(right_half)
    
    num += num_left + num_right
    stages += max(stages_left,stages_right)

    return num, stages


def group_division_testing_T2(s, threshold = 0.20):
    """
    Adapt Group Division Testing to work for T2 tests
    """

    num, stages = 0, 0
    if len(s) == 0: 
        return num, stages

    approx_num_infected = perform_T2_test(s)
    num+=1
    stages += 1

    if approx_num_infected == 0: return num, stages

    # Lesser the proportion of infected, larger the group size?
    if approx_num_infected/len(s) < threshold:
        group_size = max(1, len(s) // 32)
    else:
        group_size = max(1, len(s) // 64)
    groups = [s[i:i + group_size] for i in range(0, len(s), group_size)]
    
    # Apply conditional parallel binary splitting on each group
    for group in groups:
        if len(group) > 1:
            num_sub, stages_sub = Qtesting2_bin(group)
            num += num_sub
            stages = max(stages, stages_sub)
    
    return num, stages



def Qtesting1_bin(s, threshold = 0.2):
    """
    Aim: When infection rate is high, individual testing performs better than binary 
        and diagonal splitting. We leverage T1 tests to set a threshold which decides 
        when to switch to individual testing (optimized version).
    """
    num,stages = 0,0
    num_infected = sum(s)
    # I'm not incrementing tests/stages here bcs I redo this test at the start of both the functions I call
    # It is a matter of convenience that I choose not to pass this test's result as an argument
    
    if num_infected/len(s) > threshold and num_infected/len(s) < (1-threshold):
        temp_num, temp_stages = optimized_individual_testing(s)
        num += temp_num
        stages += temp_stages
        return num, stages

    temp_num, temp_stages = parallel_binary_splitting(s)
    num += temp_num
    stages += temp_stages

    return num,stages


def Qtesting2_bin(s, threshold = 0.2):
    """
    Aim: When infection rate is high, individual testing performs better than binary 
        and diagonal splitting. We leverage T1 tests to set a threshold which decides 
        when to switch to individual testing (optimized version).
    """
    num,stages = 0,0
    num_infected = perform_T2_test(s)
    # I'm not incrementing tests/stages here bcs I redo this test at the start of both the functions I call
    # It is a matter of convenience that I choose not to pass this test's result as an argument
    
    if num_infected/len(s) > threshold:
        temp_num, temp_stages = individual_testing(s)
        num += temp_num
        stages += temp_stages
        return num, stages

    temp_num, temp_stages = parallel_binary_splitting_T2(s)
    num += temp_num
    stages += temp_stages

    return num,stages


def Qtesting1(s):
    # Return our best performing model :)
    return group_division_testing(s)


def Qtesting2(s):
    return group_division_testing_T2(s)



def pool_communities(communities):
    pooled_communities = []
    
    # Iterate over the communities in steps of 2
    for i in range(0, len(communities), 2):
        # Pool every two consecutive communities
        pooled_community = communities[i]
        if i + 1 < len(communities):
            pooled_community = np.concatenate((communities[i], communities[i + 1]))
        pooled_communities.append(pooled_community)
    
    return np.array(pooled_communities)


def Qtesting1_comm_aware(s, communities, threshold = 1/32):
    '''
    s (np.array): binary string of infection status
    communities (list of lists): the community information
    '''
    num_tests = 0
    stages = 0

    # Initial test to determine the presence of any infection
    num_infected = np.sum(s)
    num_tests += 1
    stages += 1
    
    if num_infected == 0 or num_infected == len(s):
        return num_tests, stages
    
    # Pooling communities until a certain size
    while len(communities[0]) / len(s) <= threshold:
        communities = pool_communities(communities)

    communities = np.array(communities)
    s = np.array(s)

    for community in communities:
        temp_tests, temp_stages = Qtesting1_bin(s[community])
        num_tests += temp_tests
        stages = max(temp_stages, stages)
    
    return num_tests, stages


def Qtesting1_comm_aware_bin(s, communities):
    '''
    s(np.array): binary string of infection status
    communities(list): the community information
    '''
    num_tests = 0
    stages = 0
    communities = np.array(communities)
    s = np.array(s)

    for community in communities:
        num_tests += 1
        if sum(community) == 0: continue
        if sum(community)/len(community) <= 1/128:
            temp_tests, temp_stages = Qtesting1_bin(s[community])
            num_tests += temp_tests
        else:
            temp_tests, temp_stages = Qtesting1(s[community])
            num_tests += temp_tests
        stages = max(temp_stages,stages)

    return num_tests, stages



def Qtesting2_comm_aware(s, communities, threshold = 1/32):
    '''
    s (np.array): binary string of infection status
    communities (list of lists): the community information
    '''
    num_tests = 0
    stages = 0

    # Initial test to determine the presence of any infection
    num_infected = perform_T2_test(s)
    num_tests += 1
    stages += 1
    
    if num_infected == 0:
        return num_tests, stages
    
    # Pooling communities until a certain size
    while len(communities[0]) / len(s) <= threshold:
        communities = pool_communities(communities)

    communities = np.array(communities)
    s = np.array(s)

    for community in communities:
        temp_tests, temp_stages = Qtesting2_bin(s[community])
        num_tests += temp_tests
        stages = max(temp_stages, stages)
    
    return num_tests, stages\
    

def create_communities(N, M):
    """
    Function to create communities with approximately equal number of nodes.

    Parameters:
    N (int): Total number of nodes.
    M (int): Number of communities.

    Returns:
    List[List[int]]: 2D array of indices for each community.
    """
    communities = []
    nodes_per_community = N // M
    extra_nodes = N % M

    start = 0
    for i in range(M):
        end = start + nodes_per_community + (1 if i < extra_nodes else 0)
        communities.append(list(range(start, end)))
        start = end

    return communities