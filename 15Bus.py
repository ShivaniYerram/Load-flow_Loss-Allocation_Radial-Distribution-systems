import numpy as np
import math
import matplotlib.pyplot as plt
nBus = 14
linedata = [[1, 0, 1, 1.35309, 1.32349, 63],
                [2, 1, 2, 1.17024, 1.14464, 100],
                [3, 2, 3, 0.84111, 0.82271, 200],
                [4, 3, 4, 1.52348, 1.02760, 63],
                [5, 1, 8, 2.01317, 1.35790, 200],
                [6, 8, 9, 1.68671, 1.13770, 200],
                [7, 1, 5, 2.55727, 1.72490, 100],
                [8, 5, 6, 1.08820, 0.73400, 100],
                [9, 5, 7, 1.25143, 0.84410, 63],
                [10, 2, 10, 1.79553, 1.21110, 200],
                [11, 10, 11, 2.44845, 1.65150, 100],
                [12, 11, 12, 2.01317, 1.35790, 63],
                [13, 3, 13, 2.23081, 1.50470, 100],
                [14, 3, 14, 1.19702, 0.80740, 200]]

MVA = 100
KV = 11
Zb = (KV**2)/MVA
for i in range(len(linedata)):
    x = linedata[i][5]
    linedata[i][5] = math.cos(0.70) * x
    linedata[i].append(math.sin(0.70) * x)

power = []
for i in linedata:
    p = complex(i[5], i[6])
    power.append(p/(1000*MVA))
for i in linedata:
    i[3] = i[3]/Zb
    i[4] = i[4]/Zb

pathmatrix = []
for j in range(nBus+1):
    temp = []
    for i in linedata:
        if i[1] == j:
            temp.append(i[2])
    pathmatrix.append(temp)


def path_f(fr, to, path=[]):
    x = pathmatrix[fr]
    path = [fr]
    for i in range(len(x)):
        p=[]
        y = pathmatrix[x[i]]
        path += [path_f(x[i], to, p)]
    return path


bibc = []
for i in range(nBus+1):
    x = [0]*(nBus+1)
    bibc.append(x)

for i in range(1, nBus+1):
    irregular_list = path_f(i, nBus)

    t1 = lambda irregular_list : [element for item in irregular_list for element in t1(item)] if type(irregular_list)\
                                                                                     is list else [irregular_list]

    for j in t1(irregular_list):
        bibc[i][j] = 1

bibc1 = np.delete(bibc, 0, 0)
bibc2 = np.delete(bibc1, 0, 1)


def pathtoend(start, end, path=[]):
    path = path + [start]
    if start == end:
        return path
    for node in pathmatrix[start]:
        if node not in path:
            newpath = pathtoend(node, end, path)
            if newpath:
                return newpath
    return None


bcbv = []
for i in range(nBus+1):
    x = [0]*(nBus+1)
    bcbv.append(x)


def branchNumber(fr, to):
    for j in linedata:
        if j[1] == fr and j[2] == to:
            return [j[0], j[3], j[4]]
    return 0


for i in range(1, nBus+1):
    p = []
    t = pathtoend(0, i, p)
    for j in range(len(t)-1):
        branchNum = branchNumber(t[j], t[j+1])
        if branchNum:
            bcbv[i][branchNum[0]] = complex(branchNum[1], branchNum[2])


bcbv1 = np.delete(bcbv, 0, 0)
bcbv2 = np.delete(bcbv1, 0, 1)
tolerance = 1
dlf = np.dot(bcbv2, bibc2)
loadCurrents = [0]*nBus
voltageTemp = [0]*nBus
v1 = [complex(1, 0)] * nBus
numberIterations = 0
voltageIterations = []
while tolerance > 0.0001:
    toleranceIterations = []
    for i in range(nBus):
        loadCurrents[i] = (np.conj(power[i])) / np.conj(v1[i])
    delV = np.dot(dlf, loadCurrents)
    for j in range(nBus):
        voltageTemp[j] = complex(1, 0) - delV[j]
        toleranceIterations.append(abs(voltageTemp[j] - v1[j]))
    m = max(toleranceIterations)
    if m < 0.0001:
        tolerance = m
    for k in range(nBus):
        v1[k] = voltageTemp[k]
    numberIterations += 1
print('Voltage After Converged Load flow in ', numberIterations, 'iterations')
for i in range(len(v1)):
    print(i+2, abs(v1[i]))

for i in range(len(loadCurrents)):
    loadCurrents[i] = np.conj(power[i])/np.conj(v1[i])
branchCurrents = np.dot(bibc2, loadCurrents)
print('Load Currents After Converged Load flow')
for i in range(len(loadCurrents)):
    print(i+2, loadCurrents[i])
print(abs(min(v1)))
print(abs(max(v1)))


'''Loss Allocation Using Exact Method'''

lossAllocMatrix = []
for i in range(1, nBus+1):
    temp = []
    irregular_list = path_f(i, nBus)

    t1 = lambda irregular_list : [element for item in irregular_list for element in t1(item)] if type(irregular_list)\
                                                                                     is list else [irregular_list]
    temp.append(linedata[i-1][0])
    temp.append(linedata[i-1][1])
    temp.append(linedata[i-1][2])
    temp.append((t1(irregular_list)))
    lossAllocMatrix.append(temp)


def find_nodes(br):
    for i in linedata:
        if i[0] == br:
            return i[1],i[2]


P = 0
for i in linedata:
    x = i[0]
    ILoss = branchCurrents[x-1]
    PLoss = i[3] * (ILoss**2)
    P += abs(PLoss)
P = P*1000*MVA

powerLossAllocation = []
v2 = []
v2.append(complex(1, 0))
for i in v1:
    v2.append(i)
nodesAfterBranchMatrix = []
for o in range(nBus):
    x = [0]*nBus
    nodesAfterBranchMatrix.append(x)

for j in lossAllocMatrix:
    for k in j[3]:
        nodesAfterBranchMatrix[j[0]-1][k-1] = k


for i in range(nBus):
    voltageTerm = 0
    f = 0
    currentterm = 0
    for j in range(1, nBus+1):
        if nodesAfterBranchMatrix[j-1][i]:
            x, y = find_nodes(j)
            voltageTerm = np.conj(v2[x] - v2[y])
            currentterm = loadCurrents[nodesAfterBranchMatrix[j-1][i]-1]
            f += (voltageTerm*currentterm)

    powerLossAllocation.append(f.real*1000*MVA)

print('Loss Allocated to each node')
for i in range(len(powerLossAllocation)):
    print(i+2, powerLossAllocation[i])
print('Sum of losses allocated ', sum(powerLossAllocation))
print('Total power Loss', P)

nodes = []
absVoltages = []
absVoltages.append(1)
for i in v1:
    absVoltages.append(abs(i))
for i in range(nBus+1):
    nodes.append(i+1)
plt.ylim(0.8,1.0)
plt.ylabel('Voltage(p.u.)')
plt.xlabel('Node(Bus)')
plt.title('Load Flow Analysis')
plt.stem(nodes, absVoltages)
plt.show()

