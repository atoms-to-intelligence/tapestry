import random
import torch
import numpy as np
# import torch.nn as nn
S = 40 # number of samples
T = 16 # number of tests
d = 16  # number of samples pooled per test
k = 4  # max number of positive cases in the sample


def create_mat():
    A = np.zeros((T, S), dtype=np.uint8)

    for t in range(T):
        pool_list = []
        for i in range(d):
            r = random.randrange(0, S)
            while r in pool_list:
                r = random.randrange(0, S)
            pool_list.append(r)
            A[t, r] = 1
    return A


def increment(x, num_ones =0):
    """
        Allows enumeration of all sparse x (which has at most k ones).
        Modifies x to the next sparse x in the lexicographic ordering (with some caveats)
    """
    if num_ones == 0:
        x[0] = 1
        return x, 1

    first_one = 0
    while x[first_one] != 1:
        first_one +=1
    # print("first one", first)
    end_first_run_one = first_one
    while end_first_run_one + 1 < len(x) and x[end_first_run_one+1] != 0:
        end_first_run_one += 1
    if end_first_run_one + 1 == len(x):
        num_ones += 1
        for i in range(num_ones):
            x[i] = 1
        for i in range(num_ones, len(x)):
            x[i] = 0
        return x, num_ones
            
    run_len = end_first_run_one - first_one + 1
    x[end_first_run_one] = 0
    x[end_first_run_one+1] = 1
    for i in range(end_first_run_one + 1):
        x[i] = 0
    for i in range(run_len-1):
        x[i] = 1
    return x, num_ones

def check(A):
    """
        Check whether different sparse samples are mapped to the same test results
    """
    outs = {}
    x = [0]*S
    num_ones = 0
    while num_ones <= k:
        x_v = np.array(x, dtype=np.uint8)
        y = str(A.dot(x_v).tolist())
        # print(x_v, y)
        if y not in outs.keys():
            outs[y] = list(x)
        else:
            # print(outs[y], x_v, y)
            return False
        x, num_ones = increment(x, num_ones)
    return True

checked = False
while not checked:
    A = create_mat()
    checked = check(A)
    
    
print(A)




# i = 0
# x = [0,0,0,0,0,0]
# num_ones = 0
# print(i,x)
# while i <input_space_size - 1:
#     # print(i, x, num_ones)
#     x, num_ones = increment(x, num_ones)
#     print(i,x)
#     i+= 1




# def gen_samples(n, k, A):
#     B = np.zeros((n,S), dtype=np.uint8)
#     for i in range(n):
#         pool_list = []
#         for j in range(k):
#             r = random.randrange(0, S)
#             while r in pool_list:
#                 r = random.randrange(0, S)
#             pool_list.append(r)
#             B[i, r] = 1
#     return B, B.dot(A.transpose())


# mlp_decoder = torch.nn.Sequential(
#     nn.Linear(T, S//2),
#     nn.ReLU(),
#     nn.Linear(S//2, 2*S),
#     nn.ReLU(),
#     nn.Linear(2*S, S),
#     nn.ReLU(),
#     nn.Linear(S, S),
#     nn.ReLU(),
#     nn.Linear(S, S),
#     nn.Sigmoid()
# )


# criterion = nn.MSELoss()
# optimizer = torch.optim.SGD(mlp_decoder.parameters(), lr = 0.1)

# mlp_decoder.train()


# for i in range(10000):
#     optimizer.zero_grad()
#     y_train, x_train = gen_samples(500, k, A)
#     # print('xtrain, ytrain', x_train, y_train)
#     y_train, x_train = torch.FloatTensor(y_train), torch.FloatTensor(x_train)
#     out = mlp_decoder(x_train)
#     # print('out', out)
#     loss = criterion(y_train, out)
#     loss.backward()
#     optimizer.step()
#     out = out.round()
#     # print('out rounded', out)
#     comp = (out == y_train).float()
#     # print('comp',comp)
#     # print(comp.mean()*100)
#     comp = comp.sum(dim=1)
#     # print('summed comp',comp)
#     comp = (comp == 1.0*S).float()
#     print("acc",comp.mean()*100)
#     print(loss)


