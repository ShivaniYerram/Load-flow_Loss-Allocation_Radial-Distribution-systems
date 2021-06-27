import numpy as np
import matplotlib.pyplot as plt

linedata = [[1, 0, 1, 0.0922, 0.0470, 100, 60],
            [2, 1, 2, 0.4930, 0.2511, 90.00, 40.00],
            [3, 2, 3, 0.3660, 0.1864, 120.00, 80.00],
            [4, 3, 4, 0.3811, 0.1941, 60.00, 30.00],
            [5, 4, 5, 0.8190, 0.7070, 60, 20],
            [6, 5, 6, 0.1872, 0.6188, 200, 100],
            [7, 6, 7, 0.7114, 0.2351, 200, 100],
            [8, 7, 8, 1.0300, 0.7400, 60, 20],
            [9, 8, 9, 1.0440, 0.7400, 60, 20],
            [10, 9, 10, 0.1966, 0.0650, 45, 30],
            [11, 10, 11, 0.3744, 0.1238, 60, 35],
            [12, 11, 12, 1.4680, 1.1550, 60, 35],
            [13, 12, 13, 0.5416, 0.7129, 120, 80],
            [14, 13, 14, 0.5910, 0.5260, 60, 10],
            [15, 14, 15, 0.7463, 0.5450, 60, 20],
            [16, 15, 16, 1.2890, 1.7210, 60, 20],
            [17, 16, 17, 0.7320, 0.5740, 90, 40],
            [18, 1, 18, 0.1640, 0.1565, 90, 40],
            [19, 18, 19, 1.5042, 1.3554, 90, 40],
            [20, 19, 20, 0.4095, 0.4784, 90, 40],
            [21, 20, 21, 0.7089, 0.9373, 90, 40],
            [22, 2, 22, 0.4512, 0.3083, 90, 50],
            [23, 22, 23, 0.8980, 0.7091, 420, 200],
            [24, 23, 24, 0.8960, 0.7011, 420, 200],
            [25, 5, 25, 0.2030, 0.1034, 60, 25],
            [26, 25, 26, 0.2842, 0.1447, 60, 25],
            [27, 26, 27, 1.0590, 0.9337, 60, 20],
            [28, 21, 28, 0.8042, 0.7006, 120,70],
            [29, 28, 29, 0.5075, 0.2585, 200, 600],
            [30, 29, 30, 0.9744, 0.9630, 150, 70],
            [31, 30, 31, 0.3105, 0.3619, 210, 100],
            [32, 31, 32, 0.3410, 0.5302, 60, 40],
            ]

nBus = 32

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