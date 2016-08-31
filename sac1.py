import sys
from igraph import *
import csv
import math
import random
random.seed(0)

numOfVertices = 324
g = Graph(numOfVertices)
header = []

def create_graph():
	attrlist_file = "data/fb_caltech_small_attrlist.csv"
	edgelist_file = "data/fb_caltech_small_edgelist.txt"
	with open(attrlist_file, 'rb') as attributes_files:
		attributes = csv.reader(attributes_files)
		isHeader = True
		vertex = 0
		for row in attributes:
			if isHeader:
				for value in row:
					header.append(value)
				isHeader = False
			else:
				for index in range(0, len(header), 1):
					g.vs[vertex][header[index]] = float(row[index])
				vertex = vertex + 1
	with open(edgelist_file, 'rb') as edges_files:
		edges = csv.reader(edges_files, delimiter = " ")
		for row in edges:
			edge = (int(row[0]), int(row[1]))
			g.add_edges(edge)
			e = g.get_eid(int(row[0]), int(row[1]))
			g.es[e]["weight"] = 1.0

def getInitialCommunities():
	communities = {}
	for index in range(0, numOfVertices, 1):
		communities[index] = [index]
	return communities

def getSummationOfLinks(vertex, community):
	sum = 0
	for v in community:
		if g.are_connected(vertex, v):
			edge = g.get_eid(vertex, v)
			edgeWeight = g.es[edge]["weight"]
			sum = sum + edgeWeight
	return sum

def getSumOfEdgeWeights(vertex):
	sum = 0
	incident_edges = g.incident(vertex)
	for e in incident_edges:
		sum = sum + g.es[e]["weight"]
	return sum

def getEdgeWeightSum(community):
	sum = 0
	for v in community:
		#sum = sum + g.degree(v)
		sum = sum + getSumOfEdgeWeights(v)
	return sum

def getSimilaritySum(vertex, community):
	sum = 0
	vertex_attr = g.vs[vertex]
	for v in community:
		numerator = 0
		den_vertex = 0
		den_v = 0
		v_attr = g.vs[v]
		for h in header:
			numerator = numerator + (vertex_attr[h] * v_attr[h])
			den_vertex = den_vertex + (vertex_attr[h] * vertex_attr[h])
			den_v = den_v + (v_attr[h] * v_attr[h])
		denominator = math.sqrt(den_vertex) * math.sqrt(den_v)
		if denominator != 0:
			sum = sum + (float(numerator) / denominator)
	return sum

def sac_1_Algorithm(alpha):
	numOfEdges = len(g.get_edgelist())
	communities = getInitialCommunities()
	if numOfEdges > 0:
		maxIterations = 15
		vertices_list = list(range(0, numOfVertices, 1))
		while maxIterations > 0:
			maxIterations = maxIterations - 1
			isStabilized = True
			random.shuffle(vertices_list)
			for vertex in vertices_list:
				community_max_gain = -1
				max_gain = -1
				community_containing_vertex = -1
				for community in communities.keys():
					if vertex in communities[community]:
						communities[community].remove(vertex)
						community_containing_vertex = community
					if alpha == 0.0:
						Q_Newman = 0
					else:
						sum_links = getSummationOfLinks(vertex, communities[community])
						edge_weight_sum = getEdgeWeightSum(communities[community])
						incident_edge_weight_sum = getSumOfEdgeWeights(vertex)
						Q_Newman = (1.0 / (2 * numOfEdges)) * (sum_links - ((float(incident_edge_weight_sum) / (2 * numOfEdges)) * edge_weight_sum))
					if alpha == 1.0:
						Q_Attr = 0
					else:
						#Q_Attr = getSimilaritySum(vertex, communities[community])
						community_size = 1
						if len(communities[community]) > 0:
							community_size = len(communities[community])
						Q_Attr = getSimilaritySum(vertex, communities[community]) / community_size
					Q = (alpha * Q_Newman) + ((1 - alpha) * Q_Attr)
					if Q > 0 and Q > max_gain:
						max_gain = Q
						community_max_gain = community		
				if community_max_gain != -1:
					isStabilized = False
					communities[community_max_gain].append(vertex)
					com = community_containing_vertex
					if len(communities[com]) == 0:
						del communities[com]
				else:
					com = community_containing_vertex
					communities[com].append(vertex)
			if isStabilized:
				break
	return communities

def writeCommunitiesToFile(alpha, communities):
	if alpha == 0.0:
		name = '0'
	elif alpha == 0.5:
		name = '5'
	elif alpha == 1.0:
		name = '1'
	else:
		name = str(alpha)
	fileName = "communities_" + name +".txt"
	file = open(fileName, 'w')
	for community in communities.keys():
		length = len(communities[community])
		for node in communities[community]:
			length = length - 1;
			if(length > 0):
				file.write(str(node) + ", ")
			else:
				file.write(str(node) + "\n")
	file.close()

def main(args):
	if len(args) == 1:
		global numOfVertices
		communities = {}
		prev_communities = {}
		alpha = float(args[0])
		if alpha >= 0 and alpha <= 1:
			create_graph()
			iteration = 1
			maximumIterations = 15
			while iteration <= maximumIterations:
				if iteration != 1:
					c = 0
					membership = list(range(0, numOfVertices, 1))
					for community in communities.keys():
						for index in communities[community]:
							membership[index] = c
						c = c + 1
					numOfVertices = len(communities)
					g.contract_vertices(membership, combine_attrs = mean)
					g.simplify(combine_edges = sum)
				communities = sac_1_Algorithm(alpha)
				if len(communities) == len(prev_communities):
					break
				else:
					if iteration == 1:
						i = 0
						for com in communities.keys():
							prev_communities[i] = communities[com]
							i = i + 1
					else:
						newCom = {}
						i = 0
						for com in communities.keys():
							com_nodes = communities[com]
							mergedlist = []
							for node in com_nodes:
								mergedlist.extend(prev_communities[node])
							newCom[i] = mergedlist
							i = i + 1
						prev_communities = newCom.copy()
				iteration = iteration + 1
			writeCommunitiesToFile(alpha, prev_communities)
		else:
			print("Alpha should be between 0 and 1.")
	else:
		print("Error in number of arguments.")


if __name__ == '__main__':
	main(sys.argv[1:])