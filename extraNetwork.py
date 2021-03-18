import networkx
import re
import random
import itertools
from numpy import cumsum

transitionMat=[[1,0],[0,1]]
typeDistribution=[0.5,0.5]

def stocasticBlockModel(n,typeDistribution,transitionMat):
	types=[(random.random() < cumsum(typeDistribution)).sum()-1 for x in range(n)]
	pr1=networkx.Graph([item for item in itertools.combinations(range(n),2) if random.random()<transitionMat[types[item[0]]][types[item[1]]]])
	return pr1


def __cutHelper__(text):
	txt1=text.group().split('"')
	return '"'.join([txt1[0],txt1[0]+txt1[1],txt1[2]])

def fixNamingProblem(filename):
	g=re.sub("\d+ \"[^\n]+\"",__cutHelper__,open(filename,'r').read())
	open(filename,'w').write(g)
	return

def barabasi_albert_graph_modified(n, m,power, seed=None):
	if seed is not None:
		random.seed(seed)
	G=networkx.empty_graph(m)
	targets=list(range(m))
	repeated_nodes=[]
	source=m
	while source<n:
		G.add_edges_from(zip([source]*m,targets))
		for item in targets:
			 repeated_nodes.extend([item]*(len(G[item])**power-(len(G[item])-1)**power))
		repeated_nodes.extend([source]*m)
		targets = networkx.random_graphs._random_subset(repeated_nodes,m)
		source += 1
	return G

def snowballSample(G,seeds,hops):
	nodesToAdd=set(seeds) 
	#Not very efficient but should work for the small graphs in the exercise
	for i in range(hops):
		nodesToAdd.update([x for y in nodesToAdd for x in G[y]])
	return G.subgraph(list(nodesToAdd))

class graphEdgeDisplay(networkx.Graph):
	def __str__(self):
		return str(self.edges())
	def __repr__(self):
		return str(self.edges())

import itertools as it
def checkMotifs(G,s):
	if (s==3):
		edgeLists=[[[1,2],[1,3]],[[1,2],[1,3],[2,3]]]
	elif (s==4):
		edgeLists=[[[1,2],[1,3],[1,4]]]
		edgeLists.append([[1,2],[1,3],[1,4],[2,3]])
		edgeLists.append([[1,2],[1,3],[1,4],[2,3],[3,4]])
		edgeLists.append([[1,2],[1,3],[1,4],[2,3],[3,4],[2,4]])
	else:
		raise networkx.NetworkXNotImplemented('Size of motif must be 3 or 4')
	listOfMotifs=[graphEdgeDisplay(x) for x in edgeLists]
	counter=[0 for i in range(len(edgeLists))]
	edgeCount=[s-1,s]
	for n in G:
		for nodes in it.combinations(G[n],s-1):
			h=0
			s1=set(nodes)
			ed1=s-1+len([x for x in it.combinations(nodes,2) if x[0] in G[x[1]]])
			counter[ed1-s+1]+=1
	
	for i in range(len(listOfMotifs)):
		counter[i]=counter[i]/sum(1 for x in listOfMotifs[i] if len(listOfMotifs[i][x])==s-1)
	if (s==4):
	# 	#Need extra motifs u's and squares
		listOfMotifs.append(graphEdgeDisplay([[1,2],[2,3],[3,4]]))
		count=0
		for nod1 in G:
			for nod2 in G[nod1]:
				c3=len(set(G[nod1]).intersection(G[nod2]))
				c1=len(G[nod1])-c3-1
				c2=len(G[nod2])-c3-1				
				count+=c3*(c3-1)+(c1+c2)*c3+c1*c2
		counter.append(count/2-12*counter[3]-6*counter[2]-2*counter[1])
		listOfMotifs.append(graphEdgeDisplay([[1,2],[2,3],[3,4],[1,4]]))
		G1=G
		count=0
		for nod1 in G1:
			for nod2 in G1[nod1]:
				s2=set(G1[nod2])
				for nod3 in G1[nod1]:
					if nod2==nod3:
						continue
					count+=len(s2.intersection(G[nod3]))-1
		counter.append((count-8*counter[2]-24*counter[3])/8)
		counter[-2]-=counter[-1]*4
	# 	counts
		
	return zip(listOfMotifs,counter)
			
	
