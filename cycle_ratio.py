import time
import networkx as nx
import itertools


def read_graph(path: str):
    graph = nx.Graph()
    file = open(path)
    while 1:
        lines = file.readlines(10000)
        if not lines:
            break
        for line in lines:
            line = line[:-1]
            graph.add_edge(int(line[:6]), int(line[7:]))
    file.close()
    return graph


def node_edge_count(graph: nx.Graph):
    return graph.number_of_nodes(), graph.number_of_edges()


class CycleRatio:
    def __init__(self, graph: nx.Graph) -> None:
        self.graph = graph
        self.NodeNum, _ = node_edge_count(graph=graph)
        self.DEF_IMPOSSLEN = self.NodeNum + 1
        self.SmallestCycles = set()
        self.NodeGirth = dict()
        self.NumSmallCycles = 0
        self.CycLenDict = dict()
        self.CycleRatio = {}
        self.SmallestCyclesOfNodes = {}
        self.Coreness = nx.core_number(graph)

    def prune_graph(self):
        removeNodes = set()
        for i in self.graph.nodes():
            self.SmallestCyclesOfNodes[i] = set()
            self.CycleRatio[i] = 0
            if self.graph.degree(i) <= 1 or self.Coreness[i] <= 1:
                self.NodeGirth[i] = 0
                removeNodes.add(i)
            else:
                self.NodeGirth[i] = self.DEF_IMPOSSLEN

        # self.graph.remove_nodes_from(removeNodes)
        self.NumNode = self.graph.number_of_nodes()

        for i in range(3, self.graph.number_of_nodes()+2):
            self.CycLenDict[i] = 0

    def my_all_shortest_paths(self, source: int, target: int):
        pred = nx.predecessor(self.graph, source)
        if target not in pred:
            raise nx.NetworkXNoPath(
                f"Target {target} cannot be reached" f"from given sources"
            )
        sources = {source}
        seen = {target}
        stack = [[target, 0]]
        top = 0
        while top >= 0:
            node, i = stack[top]
            if node in sources:
                yield [p for p, n in reversed(stack[: top + 1])]
            if len(pred[node]) > i:
                stack[top][1] = i + 1
                next = pred[node][i]
                if next in seen:
                    continue
                else:
                    seen.add(next)
                top += 1
                if top == len(stack):
                    stack.append([next, 0])
                else:
                    stack[top][:] = [next, 0]
            else:
                seen.discard(node)
                top -= 1

    def getandJudgeSimpleCircle(self, objectList):
        numEdge = 0
        for eleArr in list(itertools.combinations(objectList, 2)):
            if self.graph.has_edge(eleArr[0], eleArr[1]):
                numEdge += 1
        if numEdge != len(objectList):
            return False
        else:
            return True

    def getSmallestCycles(self):
        NodeList = list(self.graph.nodes())
        NodeList.sort()
        # setp 1
        curCyc = list()
        for ix in NodeList[:-2]:  # v1
            if self.NodeGirth[ix] == 0:
                continue
            curCyc.append(ix)
            for jx in NodeList[NodeList.index(ix) + 1: -1]:  # v2
                if self.NodeGirth[jx] == 0:
                    continue
                curCyc.append(jx)
                if self.graph.has_edge(ix, jx):
                    for kx in NodeList[NodeList.index(jx) + 1:]:  # v3
                        if self.NodeGirth[kx] == 0:
                            continue
                        if self.graph.has_edge(kx, ix):
                            curCyc.append(kx)
                            if self.graph.has_edge(kx, jx):
                                self.SmallestCycles.add(tuple(curCyc))
                                for i in curCyc:
                                    self.NodeGirth[i] = 3
                            curCyc.pop()
                curCyc.pop()
            curCyc.pop()

        ResiNodeList = []  # Residual Node List
        for nod in NodeList:
            if self.NodeGirth[nod] == self.DEF_IMPOSSLEN:
                ResiNodeList.append(nod)
        if len(ResiNodeList) == 0:
            return
        else:
            visitedNodes = dict.fromkeys(ResiNodeList, set())
            for nod in ResiNodeList:
                if self.Coreness[nod] == 2 and self.NodeGirth[nod] < self.DEF_IMPOSSLEN:
                    continue
                for nei in list(self.graph.neighbors(nod)):
                    if self.Coreness[nei] == 2 and self.NodeGirth[nei] < self.DEF_IMPOSSLEN:
                        continue
                    if not nei in visitedNodes.keys() or not nod in visitedNodes[nei]:
                        visitedNodes[nod].add(nei)
                        if nei not in visitedNodes.keys():
                            visitedNodes[nei] = set([nod])
                        else:
                            visitedNodes[nei].add(nod)
                        if self.Coreness[nei] == 2 and self.NodeGirth[nei] < self.DEF_IMPOSSLEN:
                            continue
                        self.graph.remove_edge(nod, nei)
                        if nx.has_path(self.graph, nod, nei):
                            for path in self.my_all_shortest_paths(nod, nei):
                                lenPath = len(path)
                                path.sort()
                                self.SmallestCycles.add(tuple(path))
                                for i in path:
                                    if self.NodeGirth[i] > lenPath:
                                        self.NodeGirth[i] = lenPath
                        self.graph.add_edge(nod, nei)

    def find_cycle_length_frequently(self, network: str):
        self.NumSmallCycles = len(self.SmallestCycles)
        self.cycle_lenghts = []
        for cyc in self.SmallestCycles:
            lenCyc = len(cyc)
            self.cycle_lenghts.append(lenCyc)
            self.CycLenDict[lenCyc] += 1
            for nod in cyc:
                self.SmallestCyclesOfNodes[nod].add(cyc)

        import matplotlib.pyplot as plt

        # plt.subplot(211)
        # plt.title(f'Histogram lenght of smallest cycles: {network}')
        # plt.hist(self.cycle_lenghts, color='g', bins=100)
        # plt.subplot(212)
        # plt.title(f'Scalled histogram lenght of smallest cycles: {network}')
        # plt.hist(self.cycle_lenghts, color='g', bins=100)
        # plt.yscale('log')
        # plt.savefig(f'output/Histogram lenght of smallest cycles: {network}.jpg')

    def StatisticsAndCalculateIndicators(self, network: str):
        self.find_cycle_length_frequently(network=network)
        for objNode, SmaCycs in self.SmallestCyclesOfNodes.items():
            if len(SmaCycs) == 0:
                continue
            cycleNeighbors = set()
            NeiOccurTimes = {}
            for cyc in SmaCycs:
                for n in cyc:
                    if n in NeiOccurTimes.keys():
                        NeiOccurTimes[n] += 1
                    else:
                        NeiOccurTimes[n] = 1
                cycleNeighbors = cycleNeighbors.union(cyc)
            cycleNeighbors.remove(objNode)
            del NeiOccurTimes[objNode]
            sum = 0
            for nei in cycleNeighbors:
                sum += float(NeiOccurTimes[nei]) / \
                    len(self.SmallestCyclesOfNodes[nei])
            self.CycleRatio[objNode] = sum + 1

    def printAndOutput_ResultAndDistribution(self, objectList, nameString, Outpath: str):
        addrespath = Outpath + nameString + '.txt'
        Distribution = {}

        for value in objectList.values():
            if value in Distribution.keys():
                Distribution[value] += 1
            else:
                Distribution[value] = 1

        for (myk, myv) in Distribution.items():
            Distribution[myk] = myv / float(self.NodeNum)

        rankedDict_ObjectList = sorted(
            objectList.items(), key=lambda d: d[1], reverse=True)
        fileout3 = open(addrespath, 'w')
        for d in range(len(rankedDict_ObjectList)):
            fileout3.writelines("%6d %12.6f  \n" % (
                rankedDict_ObjectList[d][0], rankedDict_ObjectList[d][1]))
        fileout3.close()
        addrespath2 = Outpath + 'Distribution_' + nameString + '.txt'
        fileout2 = open(addrespath2, 'w')
        for (myk, myv) in Distribution.items():
            fileout2.writelines("%12.6f %12.6f  \n" % (myk, myv))
        fileout2.close()

    def printAndOutput_BasicCirclesDistribution(self, myCycLenDict, nameString, Outpath: str):
        Distribution = myCycLenDict
        print('\nDistribution of SmallestBasicCycles:')
        float_allBasicCircles = float(self.NumSmallCycles)
        addrespath2 = Outpath + 'Distribution_' + nameString + '.txt'
        fileout2 = open(addrespath2, 'w')
        for (myk, myv) in Distribution.items():
            if myv > 0:
                fileout2.writelines("%10d %15d  %12.6f  \n" %
                                    (myk, myv, myv/float_allBasicCircles))
                print('len:%10d,count:%10d,ratio:%12.6f' %
                      (myk, myv, myv/float_allBasicCircles))
        fileout2.close()

        List = list(self.SmallestCycles)
        rankedSBC_Set = sorted(List, key=lambda d: len(d), reverse=True)
        addrespath3 = Outpath + nameString + '_all_smallest_basic_cycles.txt'
        fileout3 = open(addrespath3, 'w')
        for cy in rankedSBC_Set:
            fileout3.writelines("%s\n" % list(cy))
        fileout3.close()
