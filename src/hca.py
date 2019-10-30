import math
import numpy as np
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

class MinHeap:
    """ Based on Max Heap implementation in Necaise, "Data structures and algorithms in Python" (Ch 13); fixed a number of bugs though... """
    def __init__(self, registry):
        self._elements = np.array([None for _ in range(registry.NR)])
        self._count = 0
        self.registry = registry
    def __len__(self):
        return self._count
    def __str__(self):
        return str([y for y in self._elements[:self._count]])
    def __repr__(self):
        return self.__str__()
    def capacity(self):
        return len(self._elements)
    def _get(self, i):
        return self.registry.heap[self._elements[i]]
    def add(self, address):
        assert self._count < self.capacity(), "Cannot add to a full heap"
        self._elements[self._count] = address
        self._count += 1
        self._siftUp(self._count - 1)
    def extract(self):
        assert self._count > 0, "Cannot extract from an empty heap"
        address = self._elements[0]
        self._count -= 1
        self._elements[0] = self._elements[self._count]
        self._siftDown(0)
        return address
    def _delete(self, i):
        assert self._count > i, "Cannot delete index" + str(i)
        self._count -= 1
        self._elements[i] = self._elements[self._count]
        self._siftDown(i)
    def _siftUp(self, i):
        if i > 0:
            parent = (i-1) // 2
            if self._get(i) < self._get(parent): # swap
                tmp = self._elements[i]
                self._elements[i] = self._elements[parent]
                self._elements[parent] = tmp
                self._siftUp(parent)
    def _siftDown(self, i):
        left = 2 * i + 1
        right = 2 * i + 2
        smallest = i
        if left < self._count and self._get(left) <= self._get(smallest):
            smallest = left
        if right < self._count and self._get(right) <= self._get(smallest):
            smallest = right
        if smallest != i: # swap
            tmp = self._elements[i]
            self._elements[i] = self._elements[smallest]
            self._elements[smallest] = tmp
            self._siftDown(smallest)

def myHash(N, idx1, idx2):
    if idx1 == idx2:
        return None
    if idx1 < idx2:
        return idx1 * (2 * N - 1) + idx2 - ((idx1 + 1) ** 2 // 2 + (idx1 + 2) // 2)
    else:
        return idx2 * (2 * N - 1) + idx1 - ((idx2 + 1) ** 2 // 2 + (idx2 + 2) // 2)

class DNode:
    """ A class for a node in a dendrogram for agglomerative clustering, assuming ultrametricity """

    def __init__(self, idx = None, children = [], dist = None):
        """ Initialise a (parent) node by the children it joins.
            Children already needs to have been initialised
            Use dist to indicate the distance between each child and this parent."""
        self.parent = None
        self.idx = idx
        self.children = children
        self.dist = dist
        if children:
            self.leaves = 0
            self.refs = []
            for child in children: # link each child to this parent
                child.parent = self
                self.leaves += child.leaves
                self.refs.extend(child.refs)
        else:
            self.leaves = 1
            self.refs = [idx]

    def nChildren(self):
        if self.children:
            return len(self.children)
        else:
            return 0

    def makeString(self, labels):
        """ Returns string with node (incl descendants) in a Newick style. """
        stubs = ['' for _ in range(self.nChildren())]
        label = dist = ''
        for i in range(self.nChildren()):
            stubs[i] = self.children[i].makeString(labels)
        if self.dist or self.dist == 0.0:
            dist = ':' + str(self.dist)
        if self.idx != None and self.idx < len(labels):
            label = labels[self.idx]
        if self.nChildren() == 0:
            return label + dist
        else:
            stubstr = '('
            for i in range(len(stubs) - 1):
                stubstr += stubs[i] + ','
            return stubstr + stubs[-1] + ')' + label + dist

    def __repr__(self):
        return '<' + str(self.idx) + '(' + str(self.leaves) + '|' + str(len(self.refs)) + '):' + str(self.dist) + ">" if self.dist else '<' + str(self.idx) + ">"

    def getDistances(self, sort = False):
        if self.children:
            d = [self.dist]
            for child in self.children:
                d.extend(child.getDistances() if not sort else sorted(child.getDistances(), reverse = True))
            return d
        else:
            return []

    def findNode(self, idx):
        if self.leaves == 1:
            return self if idx in self.refs else None
        else:
            for child in self.children:
                if idx in child.refs:
                    catch = child.findNode(idx)
                    if catch:
                        return catch
            return None

    def findParents(self, idx):
        if self.leaves == 1:
            return [self] if idx in self.refs else []
        else:
            for child in self.children:
                if idx in child.refs:
                    catch = child.findParents(idx)
                    if catch:
                        ret = [self]
                        ret.extend(catch)
                        return ret
            return []

    def getNodes(self, mindist = 0, internal = False):
        if self.leaves == 1:
            return [] if internal else [self]
        else:
            nodes = []
            if self.dist < mindist:
                return self
            for child in self.children:
                if child.leaves == 1 and not internal:
                    nodes.extend(child.getNodes(mindist))
                elif child.leaves > 1:
                    if child.dist >= mindist:
                        nodes.extend(child.getNodes(mindist))
                    else:
                        nodes.append(child)
            nodes.append(self)
            return nodes

class _RegistryIterator:
    def __init__(self, registry):
        assert len(registry._inuse) > 1, "Must have at least two indices to iterate"
        self._registry = registry
        self._n = len(registry._inuse)
        self._inuse = list(registry._inuse)
        self._curi = 0
        self._curj = 1
    def __iter__(self):
        return self
    def __next__(self):
        for i in range(self._curi, self._n - 1):
            self._curj = max(self._curi + 1, self._curj)
            for j in range(self._curj, self._n):
                address = self._registry.encode(self._inuse[i], self._inuse[j])
                y = self._registry._elements[address]
                if y:
                    if j < self._n - 1:
                        self._curj = j + 1
                    else:
                        self._curi = i + 1
                        self._curj = self._curi + 1
                    return y
        else:
            raise StopIteration

class Registry:

    def __init__(self, N):
        self.N = N                      # number of leaves
        self.NwP = 2 * N - 1            # number of leaves and parents
        self.NR = self.NwP ** 2 // 2 - self.NwP // 2  # number of records to keep
        self._elements = np.array([None for i in range(self.NR)])    # self.NR records
        self._inuse = set([i for i in range(N)]) # define what indices that are currently in use
        self.k = N

    def __iter__(self):
        return _RegistryIterator(self)

    def useIndex(self, idx):
        self._inuse.add(idx)

    def removeIndex(self, idx):
        self._inuse.remove(idx)

    def encode(self, idx1, idx2):
        if idx1 == idx2:
            return None
        if idx1 < idx2:
            return idx1 * (2 * self.N - 1) + idx2 - ((idx1+1)**2 // 2 + (idx1+2) // 2)
        else:
            return idx2 * (2 * self.N - 1) + idx1 - ((idx2+1)**2 // 2 + (idx2+2) // 2)

    def putDNode(self, dnode):
        if dnode.children[0].idx == None:
            dnode.children[0].idx = self.k
            self.k += 1
        if dnode.children[1].idx == None:
            dnode.children[1].idx = self.k
            self.k += 1
        idx1 = dnode.children[0].idx
        idx2 = dnode.children[1].idx
        address = self.encode(idx1, idx2)
        self._elements[address] = dnode

    def get(self, address):
        return self._elements[address]

    def set(self, address, value):
        self._elements[address] = value

    def getDNode(self, idx1, idx2):
        address = self.encode(idx1, idx2)
        if self._elements[address]:
            return self._elements[address]
        else:
            return None

    def decode(self, address):
        # address = idx1 * (2 * self.N - 1) + idx2 - ((idx1+1)**2 // 2 + (idx1+2) // 2)
        dnode = self._elements[address]
        if dnode:
            return [dnode.children[0].idx, dnode.children[1].idx]
        # idx1 = self._elements[address][0]
        # idx2 = address - idx1 * (2 * self.N - 1) + ((idx1+1)**2 // 2 + (idx1+2) // 2)
        return None

    def getClosest(self):
        closest = None
        try:
            for node in self:
                if closest == None or node.dist < closest.dist:
                    closest = node
        except:
            pass
        return closest

class Dendrogram:

    """ A class for a dendrogram for agglomerative clustering, assuming ultrametricity """

    def __init__(self):
        self.Labels = []                        # Labels, if available
        self.Data = []                          # Data, if available
        self.BaseTree = None                    # Root of base tree/dendrogram, if available
        self.defaultview = 'original'           # current, default view is 'original'
        self.views = {self.defaultview: None}   # root of each view
        self.means = {}                         # dict frozenset/refs : list/v

    def loadBaseTree(self, filename):
        """ Read file on Newick format.
            Returns an instance of a PhyloTree."""
        f = open(filename)
        string = ''.join(f)
        self.BaseTree = self.parseNewick(string)
        self.N = k = len(self.Labels)
        # re-indexed, viewable version
        nodes = self.BaseTree.getNodes()
        for i in range(len(nodes)):
            if not nodes[i].leaves == 1:  # not a child
                nodes[i].idx = k
                k += 1
        self.createView(self.defaultview, self.getNodes())

    def parseNewick(self, string):
        """ Main method for parsing a Newick string into a (phylogenetic) tree.
            Handles labels (on both leaves and internal nodes), and includes distances (if provided).
            Returns an instance of a PhyloTree. """
        if string.find(';') != -1:
            string = string[:string.find(';')]
        return self.parseNewickNode(string)

    def saveBaseTree(self, filename):
        assert self.BaseTree and self.Labels != None, "Invalid tree: cannot save"
        with open(filename, 'w') as fh:
            print(self.BaseTree.makeString(self.Labels), end="", file=fh)

    def loadData(self, filename, headers = False):
        reindex = False
        if self.Labels and self.BaseTree: # labels are already available for base tree, so data need to be re-indexed
            reindex = True
            self.Data = [None for _ in range(len(self.Labels))]
        with open(filename) as csvfile:
            r = csv.reader(csvfile, delimiter='\t')
            headers = [] if headers else None
            for row in r:
                if headers == []:
                    headers = row[1:]
                elif reindex:
                    self.Data[self.indexLabel(row[0])] = [float(y) for y in row[1:]]
                else:
                    self.Data.append([float(y) for y in row[1:]])
                    self.Labels.append(row[0])

    def setData(self, data, labels):
        assert len(data) == len(labels), "Number of data rows and labels must be the same"
        reindex = False
        if self.Labels and self.BaseTree: # labels are already available for base tree, so data need to be re-indexed
            reindex = True
            self.Data = [None for _ in range(len(self.Labels))]
        for i in range(len(data)):
            if reindex:
                idx = self.indexLabel(labels[i])
                self.Data[idx] = data[idx]
            else:
                self.Data.append(data[i])
                self.Labels.append(labels[i])

    def createBaseTree(self, d, labels = None):
        " d is the distance matrix, labels is the list of labels "
        self.N = len(d)
        self.Labels = labels or ['X' + str(i + 1) for i in range(self.N)]
        registry = Registry(self.N)
        " For each node-pair, assign the distance between them. "
        dnodes = [DNode(idx = i) for i in range(self.N)]
        for i in range(self.N):
            registry.useIndex(i)
            for j in range(i + 1, self.N):
                registry.putDNode(DNode(children = [dnodes[i], dnodes[j]], dist = d[i, j]))
        # Treat each node as a cluster, until there is only one cluster left, find the closest pair
        # of clusters, and merge that pair into a new cluster (to replace the two that merged).
        # In each case, the new cluster is represented by the node that is formed.
        closest = registry.getClosest()
        while closest:
            # So we know the closest, now we need to merge...
            x = closest.children[0]
            y = closest.children[1]
            z = closest  # use this node for new cluster z
            Nx = x.leaves  # find number of sequences in x
            Ny = y.leaves  # find number of sequences in y
            for widx in registry._inuse:  # for each node w ...
                if widx != x.idx and widx != y.idx:
                    # we will merge x and y into a new cluster z, so need to consider w (which is not x or y)
                    dxw_address = registry.encode(x.idx, widx)
                    dyw_address = registry.encode(y.idx, widx)
                    dxw = registry.get(dxw_address)  # retrieve and remove distance from D: x to w
                    dyw = registry.get(dyw_address)  # retrieve and remove distance from D: y to w
                    w = dxw.children[0] if dxw.children[0].idx == widx else dxw.children[1]
                    dzw = DNode(children = [z, w], dist = (Nx * dxw.dist + Ny * dyw.dist) / (Nx + Ny))  # distance: z to w
                    registry.set(dxw_address, None)
                    registry.set(dyw_address, None)
                    registry.putDNode(dzw)
            registry.useIndex(z.idx)
            # remove x and y from registry
            registry.removeIndex(x.idx)
            registry.removeIndex(y.idx)
            closest = registry.getClosest()
        self.BaseTree = z
        # re-indexed, viewable version
        self.createView(self.defaultview, self.getNodes())

    def getRoot(self, view = None):
        if not view:
            view = self.defaultview
        return self.views[view]

    def getNodes(self, mindist = 0, internalOnly = False, view = None):
        if not view:
            return self.BaseTree.getNodes(mindist, internalOnly)
        else:
            return self.views[view].getNodes(mindist, internalOnly)

    def getIndex(self, label):
        for i in range(len(self.Labels)):
            if self.Labels[i] == label:
                return i
        return -1

    def getNodeByName(self, label, view = None):
        idx = self.getIndex(label)
        if idx != -1:
            return self.getRoot(view).findNode(idx)
        else:
            return None

    def getParentsByName(self, label, view = None):
        idx = self.getIndex(label)
        if idx != -1:
            myroot = self.getRoot(view)
            return myroot.findParents(idx)
        else:
            return None

    def _intersect(self, L1, L2):
        return [y for y in L1 if y in L2]

    def getParentsByMembers(self, members = [], view = None):
        """ Retrieve the smallest sub-tree that contains all the members. """
        prevtrace = None
        for member in members:
            trace = self.getParentsByName(member, view)
            if prevtrace:
                prevtrace = self._intersect(prevtrace, trace)
            else:
                prevtrace = trace
        return prevtrace

    def getNodesByMembers(self, members = [], view = None):
        shared = self.getParentsByMembers(members, view)
        if shared:
            return shared[-1].getNodes()
        else:
            return None

    def setView(self, name):
        self.defaultview = name

    def createView(self, name, nodes):
        newnodes = []
        old2new = {}
        leafcnt = 0
        parcnt = len(nodes) // 2 + 1
        for i in range(len(nodes)):
            oldnode = nodes[i]
            if oldnode.leaves == 1: # leaf
                newnode = DNode(idx = leafcnt)
                leafcnt += 1
            else: # children, in the old node definitely, in the new node maybe
                # either both are nodes or none is
                if oldnode.children[0].idx in old2new and oldnode.children[1].idx in old2new:
                    newnode = DNode(idx = parcnt, children = [newnodes[old2new[oldnode.children[0].idx]], newnodes[old2new[oldnode.children[1].idx]]], dist = oldnode.dist)
                    parcnt += 1
                else:
                    newnode = DNode(idx = leafcnt)
                    leafcnt += 1
            old2new[oldnode.idx] = i
            newnode.refs = oldnode.refs
            newnodes.append(newnode)
        self.views[name] = newnode
        self.setView(name)

    def getLinkage(self, nodes):
        links = []
        for n in nodes:
            if n.leaves > 1:
                links.append([n.children[0].idx, n.children[1].idx, n.dist, n.leaves])
        return links

    def getRefs(self, nodes):
        groups = []
        for n in nodes:
            if n.leaves == 1:
                groups.append(n.refs)
        return groups

    def getMeanRefs(self, refs):
        if self.Data:
            key = frozenset(refs)
            if key in self.means:
                return self.means[key]
            m = None
            for r in refs:
                if not m:
                    m = [y for y in self.Data[r]]
                else:
                    for i in range(len(m)):
                        m[i] += self.Data[r][i]
            m = [y/len(refs) for y in m]
            self.means[key] = m
            return m
        else:
            return None

    def getLabelGroups(self, nodes):
        groups = []
        for n in nodes:
            if n.leaves == 1:
                plabels = []
                for i in range(len(n.refs)):
                    plabels.append(self.Labels[n.refs[i]])
                groups.append(plabels)
        return groups

    def getString(self, v, s, thresholds):
        assert len(s) - 1 == len(thresholds), "characters in string s must agree with thresholds + 1"
        ret = ''
        for y in v:
            for i in range(len(thresholds)):
                if y < thresholds[i]:
                    ret += s[i]
                    break
            if y >= thresholds[-1]:
                ret += s[-1]
        return ret

    def draw(self, view = None, maxlabel = 15, mode = 'label'):
        N = self.getRoot(view).getNodes()
        labs = []
        if mode == 'label':
            groups = self.getLabelGroups(N)
            for group in groups:
                label = '' if len(group) == 1 else str(len(group)) + ':'
                for i in range(len(group)):
                    label += group[i]
                    if len(label) > maxlabel:
                        label = label[0:maxlabel-2] + '..'
                labs.append(label)
        elif mode.startswith(' '):
            mvecs = []
            tmin =  99999999
            tmax = -99999999
            groups = self.getRefs(N)
            for group in groups:
                meanv = self.getMeanRefs(group)
                mvecs.append(meanv)
                for y in meanv:
                    tmin = min(tmin, y)
                    tmax = max(tmax, y)
            for meanv in mvecs:
                labs.append(self.getString(meanv, mode, [(tmax - tmin) / 4 + tmin, (tmax - tmin) / 3 + tmin, (tmax - tmin) / 2 + tmin]))
        L = self.getLinkage(N)
        Z = np.array(L)
        plt.rcParams["font.family"] = "Andale Mono"
        hierarchy.dendrogram(Z, labels=labs, orientation = 'left')

    """ ----------------------------------------------------------------------------------------
        Methods for processing files of trees on the Newick format
        ----------------------------------------------------------------------------------------"""

    def _findComma(self, string, level=0):
        """ Find all commas at specified level of embedding """
        mylevel = 0
        commas = []
        for i in range(len(string)):
            if string[i] == '(':
                mylevel += 1
            elif string[i] == ')':
                mylevel -= 1
            elif string[i] == ',' and mylevel == level:
                commas.append(i)
        return commas

    def indexLabel(self, label):
        if self.Labels == None:
            self.Labels = []
        for i in range(len(self.Labels)):
            if self.Labels[i] == label:
                return i
        self.Labels.append(label)
        return len(self.Labels) - 1

    def parseNewickNode(self, string):
        """ Utility function that recursively parses embedded string using Newick format. """
        first = string.find('(')
        last = string[::-1].find(')')  # look from the back
        if first == -1 and last == -1:  # we are at leaf
            y = string.split(':')
            node = DNode(idx = self.indexLabel(y[0]))
            if len(y) >= 2:
                node.dist = float(y[1])
            return node
        elif first >= 0 and last >= 0:
            # remove parentheses
            last = len(string) - last - 1  # correct index to refer from start instead of end of string
            embed = string[first + 1:last]
            tail = string[last + 1:]
            # find where corresp comma is
            commas = self._findComma(embed)
            if len(commas) < 1:
                raise RuntimeError('Invalid format: invalid placement of "," in sub-string "' + embed + '"')
            prev_comma = 0
            child_tokens = []
            for comma in commas:
                child_tokens.append(embed[prev_comma:comma].strip())
                prev_comma = comma + 1
            child_tokens.append(embed[prev_comma:].strip())
            y = tail.split(':')
            node = DNode()
            if len(y) >= 2:
                node.dist = float(y[1])
            node.children = []
            node.leaves = 0
            node.refs = []
            for tok in child_tokens:
                child = self.parseNewickNode(tok)
                child.parent = node
                node.leaves += child.leaves
                node.refs.extend(child.refs)
                node.children.append(child)
            return node
        else:
            raise RuntimeError('Invalid format: unbalanced parentheses in sub-string "' + string + '"')

"""

"""

def eucdist(v1, v2):
    diff = 0
    for i in range(len(v1)):
        diff += (v1[i] - v2[i])**2
    return math.sqrt(diff)

def cosdist(v1, v2):
    if len(v1) != len(v2):
        return None
    sum0 = 0
    sum1 = 0
    sum2 = 0
    for i in range(len(v1)):
        sum0 += v1[i]*v2[i]
        sum1 += v1[i]*v1[i]
        sum2 += v2[i]*v2[i]
    return 1 - (sum0 / (math.sqrt(sum1*sum2)))

def normalise(v):
    y = sum(v)
    return [vi/y for vi in v]

def euc(data):
    d = np.zeros((len(data), len(data)))
    for i in range(len(data)):
        for j in range(i + 1, len(data)):
            d[i, j] = eucdist(data[i], data[j])
    return d

def norm(data):
    ndata = [[] for _ in data]
    for i in range(len(data)):
        ndata[i] = normalise(data[i])
    return ndata

import csv

if __name__ == '__main__':
    dgram = Dendrogram()
    dgram.loadData('/Users/mikael/simhome/monocyte/test100x10.tsv', headers=True)
    d = euc(dgram.Data)
    dgram.createBaseTree(d, dgram.Labels)
    dgram.createView('focus slamf7', dgram.getNodesByMembers(['p3@SLAMF7', 'p2@NFKB2'], view = 'original'))
    dgram.createView('focus sorl1', dgram.getNodesByMembers(['p2@SORL1', 'p1@FAU'], view = 'original'))
    dgram.createView('focus ddx6', dgram.getNodesByMembers(['p3@DDX6', 'p2@CARS'], view = 'original'))
    dgram.createView('focus sorl1 dist', dgram.getNodes(mindist = 0.04, view = 'focus sorl1'))
    plt.figure(num=None, figsize=(5, 15), dpi=200, facecolor='w', edgecolor='k')
    #dgram.draw('focus ddx6', mode = 'label')
    dgram.draw('focus ddx6', mode = ' .o0')
    plt.show()

if __name__ == '__main__1':
    dgram = Dendrogram()
    dgram.loadData('/Users/mikael/simhome/monocyte/test10x10.tsv', headers=True)
    d = euc(norm(dgram.Data))
    dgram.createBaseTree(d, dgram.Labels)
    plt.figure(num=None, figsize=(5, 15), dpi=200, facecolor='w', edgecolor='k')
    dgram.draw()
    plt.show()
    dgram.saveBaseTree('/Users/mikael/simhome/monocyte/test10x10.dgram')

if __name__ == '__main__2':
    dgram = Dendrogram()
    dgram.loadData('/Users/mikael/simhome/monocyte/test10x10.tsv', headers=True)
    dgram.loadBaseTree('/Users/mikael/simhome/monocyte/test10x10.dgram')
    plt.figure(num=None, figsize=(5, 15), dpi=200, facecolor='w', edgecolor='k')
    dgram.draw()
    plt.show()

if __name__ == '__main__3':
    data = [[1.0, 1.0, 1.0],
            [2.0, 2.0, 1.9],
            [2.1, 2.1, 2.1],
            [3.1, 3.2, 3.3],
            [3.2, 3.4, 3.6],
            [3.5, 3.4, 3.3],
            [0.9, 0.9, 1.0]]
    labels = ['S1', 'S2a', 'S2b', 'S3a', 'S3b', 'S3c', 'S0']
    dgram = Dendrogram()
    dgram.createBaseTree(euc(data), labels)
    dgram.createView("small", dgram.getNodes(mindist=0.2))
    dgram.saveBaseTree('/Users/mikael/simhome/monocyte/small.dgram')
    dgram.draw()
    plt.show()

if __name__ == '__main__4':
    data = [[1.0, 1.0, 1.0],
            [2.0, 2.0, 1.9],
            [2.1, 2.1, 2.1],
            [3.1, 3.2, 3.3],
            [3.2, 3.4, 3.6],
            [3.5, 3.4, 3.3],
            [0.9, 0.9, 1.0]]
    labels = ['S1', 'S2a', 'S2b', 'S3a', 'S3b', 'S3c', 'S0']
    dgram = Dendrogram()
    dgram.loadBaseTree('/Users/mikael/simhome/monocyte/small.dgram')
    dgram.setData(data, labels)
    dgram.draw()
    plt.show()
