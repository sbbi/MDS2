#coding:utf8

"""
Author          :   Tian Gao
CreationDate    :   
Description     :
"""

class Node:
    """
    node structure in a graph
    """
    def __init__(self, title, data=None):
        self.nodeTitle = title
        self.nodeData = data
#end_class

class Edge:
    """
    edge structure in a graph
    """
    def __init__(self, startNode, endNode, title=None, weight=1, data=None):
        self.edgeTitle = title
        self.startNodeTitle = startNode.nodeTitle
        self.endNodeTitle = endNode.nodeTitle
        self.edgeWeight = weight
        self.edgeData = data
#end_class

class Graph:
    """
    adjacentNodeDict: adjacent nodes of each node
    (directed)  key: node1, value: [(node2, isNode1ToNode2)]
    (undirected)key: node1, value: [(node2, None)]

    adjacentEdgeDict: adjacent nodes of each node
    (directed)  key: node1, value: [(edge1, isIn)]
    (undirected)key: node1, value: [(edge1, None)]

    """
    def __init__(self, isDirected=False):
        self.isDirected = isDirected
        self.nodeCnt = 0
        self.edgeCnt = 0
        self.nodeDict = {}
        self.edgeDict = {}
        self.adjacentNodeDict = {}
        self.adjacentEdgeDict = {}
    #end_init

    def addNode(self, node):
        self.nodeCnt += 1
        self.nodeDict.setdefault(node.nodeTitle, node)
    #end_func

    def addEdge(self, edge):
        # if edge title is not given, the default value is the total number of edge
        edge.edgeTitle = str(self.edgeCnt) if edge.edgeTitle is None else edge.edgeTitle
        self.edgeDict.setdefault(edge.edgeTitle, edge)
        self.edgeCnt += 1

        # update adjacentNodeDict
        self.adjacentNodeDict.setdefault(edge.startNodeTitle, [])
        self.adjacentNodeDict.setdefault(edge.endNodeTitle, [])
        curStartNodeInfo = (edge.endNodeTitle, True) if self.isDirected else (edge.endNodeTitle, None)
        curEndNodeInfo = (edge.startNodeTitle, False) if self.isDirected else (edge.startNodeTitle, None)
        self.adjacentNodeDict[edge.startNodeTitle].append(curStartNodeInfo)
        self.adjacentNodeDict[edge.endNodeTitle].append(curEndNodeInfo)

        # update adjacentEdgeDict
        self.adjacentEdgeDict.setdefault(edge.startNodeTitle, [])
        self.adjacentEdgeDict.setdefault(edge.endNodeTitle, [])
        curEdgeInfo = (edge.edgeTitle, False) if self.isDirected else (edge.edgeTitle, None)
        self.adjacentEdgeDict[edge.startNodeTitle].append(curEdgeInfo)
        self.adjacentEdgeDict[edge.endNodeTitle].append(curEdgeInfo)
    #end_func
#end_class

def test():
    g = Graph(True)
    n1 = Node('n1')
    n2 = Node('n2')
    e1 = Edge(n1, n2)

    g.addNode(n1)
    g.addNode(n2)
    g.addEdge(e1)
    print g

#end_func

def main():
    test()
#end_main

if __name__ == "__main__":
   main()
#end_if


