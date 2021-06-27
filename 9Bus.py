import numpy as np
import matplotlib.pyplot as plt
nBus = 8
linedata = [[1, 0, 1, 1.632, 1.1019, 162, 96],
            [2, 1, 2, 1.088, 0.7346, 150, 138],
            [3, 2, 3, 0.5440, 0.3673, 12, 7.2],
            [4, 3, 4, 0.2720, 0.1836, 30, 18],
            [5, 1, 5, 0.5440, 0.3673, 45.6, 33.6],
            [6, 2, 6, 1.3760, 0.3896, 12, 6],
            [7, 3, 7, 2.7520, 0.7792, 180, 140.4],
            [8, 7, 8, 4.1280, 1.1688, 156, 120]]
MVA = 100
KV = 11
Zb = (KV**2)/MVA
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
        p = []
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
    t = pathtoend(0, i, p) #[1-2-3-4-7-8-9]
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
v1 = [complex(1, 0)] * nBus #1pu 0 angle
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


for i in nodesAfterBranchMatrix:
    print(i)

for i in range(nBus):
    voltageTerm = 0
    f = 0
    currentterm = 0
    for j in range(1, nBus+1):
        if nodesAfterBranchMatrix[j-1][i]:
            x, y = find_nodes(j)
            voltageTerm = np.conj(v2[x] - v2[y])
            currentterm = loadCurrents[i]
            f += (voltageTerm*currentterm)

    powerLossAllocation.append(f.real*1000*MVA)
print('Loss Allocated to each node')
for i in range(len(powerLossAllocation)):
    print(i+2,powerLossAllocation[i])
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









































