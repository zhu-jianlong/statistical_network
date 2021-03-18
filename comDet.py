import networkx as nx
import random as rd
from decimal import Decimal
from collections import Counter
def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

def hypergeoCDF(N,K,n,k):
	res=0
	for i in range(k,N):
		res+=hypergeometricDist(N,K,n,i)
	return res

def hypergeometricDist(N,K,n,k):
	n1=choose(K,k)
	n2=choose(N-K,n-k)
	n3=choose(N,n)
	return float(n1)*float(n2)/float(n3)

def testComEnrichment(com,data):
	com1=[[x.strip() for x in com if com[x]==y] for y in range(1+max(com.values()))]
	hometown_location=Counter(data[item]['hometown_location']['city'] if data[item]['hometown_location']!=None else None  for item in data)
	hometown_country=Counter(data[item]['hometown_location']['country'] if data[item]['hometown_location']!=None else None  for item in data)
	current_city=Counter(data[item]['current_location']['city'] if data[item]['current_location']!=None else None  for item in data)
	current_country=Counter(data[item]['current_location']['country'] if data[item]['current_location']!=None else None  for item in data)
	affiliations=Counter(aff['name'] for item in data for aff in data[item]['affiliations'])
	res=[]
	for i in range(len(com1)):
		hometown_location1=Counter(data[item]['hometown_location']['city'] if data[item]['hometown_location']!=None else None  for item in com1[i] if item in data)
		hometown_country1=Counter(data[item]['hometown_location']['country'] if data[item]['hometown_location']!=None else None  for item in com1[i] if item in data)
		current_city1=Counter(data[item]['current_location']['city'] if data[item]['current_location']!=None else None  for item in com1[i] if item in data)
		current_country1=Counter(data[item]['current_location']['country'] if data[item]['current_location']!=None else None  for item in com1[i] if item in data)
		affiliations1=Counter(aff['name'] for item in com1[i] if item in data for aff in data[item]['affiliations'])
		res1=[]
		if i==0:
			assert(0==1)
		for item in hometown_location1:	
			pValue=hypergeoCDF(len(data),hometown_location[item],len(com1[i]),hometown_location1[item])
			if pValue<0.001:
				res1.append(['homeTownCity',item,pValue])
		for item in hometown_country1:	
			pValue=hypergeoCDF(len(data),hometown_country[item],len(com1[i]),hometown_country1[item])
			if pValue<0.001:
				res1.append(['homeTownCountry',item,pValue])
		for item in current_city1:	
			pValue=hypergeoCDF(len(data),current_city[item],len(com1[i]),current_city1[item])
			if pValue<0.001:
				res1.append(['currentCity',item,pValue])
		for item in current_country1:	
			pValue=hypergeoCDF(len(data),current_country[item],len(com1[i]),current_country1[item])
			if pValue<0.001:
				res1.append(['currentCountry',item,pValue])
		for item in affiliations1:	
			pValue=hypergeoCDF(len(data),affiliations[item],len(com1[i]),affiliations1[item])
			if pValue<0.001:
				res1.append(['affiliations',item,pValue])
		
		
		res.append([x for x in res1 if x[1]!=None])
	return res
		
		
		
	return hometown_location



def loadInfo(filename):
	f=open(filename).read()
	data='['+'['.join(f.split('[')[1:])
	data=data.replace('null',"None")
	data=data.replace('false',"False")
	data=data.replace('true',"True")
	data1=eval(data)
	return {x['name'].strip('"'):x for x in data1}
			
def loadNetwork(filename):
    G=nx.Graph()
    f=open(filename).readlines()
    process=0
    for item in f:
        if item.startswith('uid1'):
            process=1
        elif process==1:
            G.add_edge(*item.strip().replace('"','').split(','))
    return G

def applyNGnullModel(B,param=1):
    k=B.sum(1)
    B-=param*k*k.T/k.sum()

# def applyNGnullModel(B,param=1):
#     for i in range(len(B)):
#         B[i][i]=0
#     k=[sum(x) for x in B]
#     twom=float(sum(k))
#     for i in range(len(B)):
#         for j in range(len(B)):
#             B[i][j]-=param*float(k[i]*k[j])/twom
#     return B

def applyUniformNullModel(B):
    for i in range(len(B)):
        B[i][i]=0
    p=sum(sum(x for x in B))/float(len(B)*len(B))
    for i in range(len(B)):
        for j in range(len(B)):
            B[i][j]-=p

def __modularity__(G,com1):
	deg=G.degree()
	twom=float(sum(list(dict(deg).values())))
	mod=0
	for item in com1:
		for i in item:
			for j in item:
				if i!=j:
					if G.has_edge(i,j):
						mod+=1
					mod-=float(deg[i]*deg[j])/twom
	return mod/twom
		
def modularity(G,com1):
	com2=defaultdict(list)
	for item in com1:
		com2[com1[item]].append(item)
	return __modularity__(G,com2.values())
	
from collections import defaultdict
def betweenComDetect(G1):
	G=G1.copy()
	curEdge=nx.edge_betweenness(G)
	res=[]
	mod1=-1000000
	while G.number_of_edges()>0:
		edge=Counter(curEdge).most_common()[0][0]
		G.remove_edge(*edge)
		curEdge=nx.edge_betweenness(G)
		com1=list(nx.connected_components(G))
		mod2=__modularity__(G1,com1)
		if mod2>mod1:
			mod1=mod2
			com2=com1
	return {x:i for i in range(len(com2)) for x in com2[i]}
		
    
import itertools as it

def louvainAlg(G,nullFunction=applyNGnullModel,param=None):
    nodes=list(G)
    nodeMap={nodes[i]:i for i in range(len(nodes))}
    B=nx.to_numpy_matrix(G)
    if param==None:
	    applyNGnullModel(B)
    else:
        applyNGnullModel(B,param)
    G1=nx.Graph()
    G1.add_edges_from((nodeMap[edge[0]],nodeMap[edge[1]]) for edge in G.edges())
    coms=genlouvain(B,G1)
    return {nodes[i]:coms[i] for i in range(len(coms))}
    
 

#
import numpy as np
# import weave

def genlouvain(B,G,limit=10000):
    uniY=range(len(B))
    Sb=[0 for i in range(len(B))]
    dtot=0
    M=B
    tree=[]
    h0=1
    # expr='Mnew[i,j]=M[ix_(y==uniY[i],y==uniY[j])].sum()'
    while (h0==1):
        h0=0
        # print len(M)
        for i in range(len(M)):
            M[i,i]=0
        y=np.arange(len(uniY))
        dstep=1
        index1=np.arange(len(y)) #range(len(y))
        h=1
        while (h==1):
            h=0
            # print 'i looped1 '+str(len(set(y)))
            rd.shuffle(index1)
            for item in [x for x in index1 if x in G]:
                # print 'i looped 2, '+str(len(set(y)))
                possLabels=list(set(y[i] for i in G[item]))+[y[item],]
                temp1=[M[item,y==i].sum() for i in possLabels]
                temp2=possLabels[temp1.index(max(temp1))]
                if y[item]!=temp2:
                    h=1
                    h0=1
                    y[item]=temp2
        uniY=list(set(y))
        mapping={uniY[i]:i for i in range(len(uniY))}
        tree.append({i:mapping[y[i]] for i in range(len(y))})
        # Mnew1=[[0 for i in range(len(uniY))] for j in range(len(uniY))]
        Mnew=np.zeros([len(uniY),len(uniY)])
        for i in range(len(uniY)):
            x1=(y==uniY[i])
            for j in range(i+1,len(uniY)):
                Mnew[i,j]=M[np.ix_(y==uniY[i],y==uniY[j])].sum()
                #weave.blitz(expr)
        Mnew=Mnew+Mnew.T
        # for i in range(len(y)):
        #     for j in range(i+1,len(y)):
        #         Mnew1[mapping[y[i]]][mapping[y[j]]]+=M[i,j]
        #         Mnew1[mapping[y[j]]][mapping[y[i]]]+=M[i,j]
        # assert(0==1)
        M=Mnew
        G=nx.Graph()
        # G.add_nodes_from(range(len(M)))
        for i in range(len(M)):
            for j in range(i+1,len(M)):
                if M[i,j]>0:
                    G.add_edge(i,j)
        

    #Need to unbundle mapping
    community=[0 for i in range(len(B))]
    for i in range(len(B)):
        cur=i
        for item in tree:
            cur=item[cur]
        community[i]=cur
    return community

	 
	 
	 
import numpy as np
def plotComs(coms):
	mat1=np.zeros([len(coms[0]),len(coms)])
	nodes=list(coms[0].keys())
	for i in range(len(coms)):
		for j in range(len(nodes)):
			mat1[j,i]=coms[i][nodes[j]]
	ix=np.lexsort([mat1[:,-i] for i in range(mat1.shape[1]+1)])
	mat1=mat1[ix]
	from pylab import ion,imshow
	ion()
	imshow(mat1,interpolation='nearest',aspect='auto')
	return

