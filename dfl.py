#Device-free localization library - The functions in this file are used
#to locate targets using RSS measurements from wireless networks.
#Author: Joey Wilson (joey.wilson@utah.edu)

from pylab import find, zeros, ones, array, arange, rand, polyfit,\
        mean, log10, exp, hist, argmax
from numpy.linalg import norm


#Skew-Laplace
def skewLaplace(x,a,b,mu):
    """Return the value of the Skew-Laplace PDF"""
    if x < mu:
        return a*b/(a+b)*exp(-a*(mu-x))
    else:
        return a*b/(a+b)*exp(-b*(x-mu))

def skewLaplaceCDF(x,a,b,mu):
    """Return the value of the Skew-Laplace PDF"""
    if x < mu:
        return b/(a+b)*exp(-a*(mu-x))
    else:
        return (b-a*exp(-b*(x-mu))+a)/(a+b)

def pooyanML(data,nbins=40):
    (fi,xi,p) = hist(data,nbins)
    maxInd = argmax(fi)
    mu = (xi[maxInd]+xi[maxInd+1])/2

    dataPlus = data[find(data>mu)] - mu
    dataMinus = -data[find(data<mu)] + mu

    alpha = 1.0/mean(dataMinus)
    beta = 1.0/mean(dataPlus)

    return mu,alpha,beta

def getPathLossParams(dVec,dataVec):
    """Estimate path loss params from data"""
  
    indx = dVec!=0.0
    dVecDb = 10*log10(dVec[indx])
    dataVecI = dataVec[indx]
    p = polyfit(dVecDb,dataVecI,1)
    return p

def nodeDistances(nodes):
    """Return a tuple of the distances from each node to the others."""
    N = len(nodes)
    distances = zeros(N**2)
    for n1 in xrange(N):
        for n2 in xrange(N):
            distances[n1*N+n2] = norm(nodes[n1] - nodes[n2])
    return distances

def fadeLevel(d,rss,params):
    """Determine the fade level of a particular link"""
    
    predicted = params[1]+params[0]*10*log10(d)
    return rss-predicted

def linkIndexAsym(node1,node2, numNodes):
    """Get the link index from two node indices. Assume one-way asymmetric links."""
    
    return node1*numNodes + node2

def nodeIndexAsym(linkIndex, numNodes):
    """Return the node indices for a given link index and number of nodes."""
    
    for n1 in range(0,numNodes):
        for n2 in range(0,numNodes):
            linkNum = n1*numNodes + n2
            if linkNum == linkIndex:
                return n1,n2

def linkIndexSym(node1,node2):
    """Get the link index from two node indices. Assume two-way symmetric links.

    This function converts the node indices into the link index.
    Each non-diagonal element from the vector is chosen row first, then column
    
    [-        ]
    [0 -      ]
    [1 2 -    ]
    [3 4 5 -  ]
    [6 7 8 9 -] etc...
    """

    if node1==node2:
        raise Exception, "Input parameters must not be equal"
    elif node1<node2:
        lowNode, highNode = node1, node2
    else:
        lowNode, highNode = node2, node1

    return sum(range(0,highNode)) + lowNode
    #return sum(range(1,highNode-1)) + lowNode

def nodeIndexSym(linkIndex,N):
    """Get the node indices from a link index.
    
    N is the number of total nodes in the network. Returned nodes are indexed from 0.
    """  
    counter = -1;
    for k in range(0,N):
        for p in range(0,k):
            if p!=k:
                counter+=1;
                if counter==linkIndex:
                    return [k,p]
            
def knownPath(numPoints,pathDict,letterString):
    """Returns the x,y coordinates of a target moving at constant velocity
    
    BEWARE: This only supports non-diagonal motion moving parallel to the x or y axis."""
    
    pathVecX = []
    pathVecY = []
    
    #Figure out the step size given total number of samples
    totalDist = 0
    for k in range(0,len(letterString)-1):
        currentDist = norm(array(pathDict[letterString[k]]) - array(pathDict[letterString[k+1]]))
        totalDist += currentDist
    stepSize = float(totalDist)/float(numPoints)
    
    for k in range(0,len(letterString)-1):
        point1 = pathDict[letterString[k]]
        point2 = pathDict[letterString[k+1]]
        
        #Is the object moving on x or y directions.  Diagonal not supported.
        if (point1[0] - point2[0] == 0):
            #Y direction
            if (point1[1] - point2[1] < 0):
                yline = arange(point1[1],point2[1],stepSize)
            else:
                yline = arange(point1[1],point2[1],-stepSize)
            xline = point1[0]*ones(len(yline))    
        else:
            #X direction
            if (point1[0] - point2[0] < 0):
                xline = arange(point1[0],point2[0],stepSize)
            else:
                xline = arange(point1[0],point2[0],-stepSize)
            yline = point1[1]*ones(len(xline))
        
        for k in range(0,len(xline)):
            pathVecX.append(xline[k])
            pathVecY.append(yline[k])
        
    return array([pathVecX[:numPoints],pathVecY[:numPoints]]).T

