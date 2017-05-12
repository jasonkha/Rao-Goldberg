import numpy as np
from scipy.linalg import eig
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans
import math
import pylab as pl
def generate_laplacian(vertices):
    laplacian=[]
    for i in range(len(vertices)):
        curr=vertices[i]
        norms=np.linalg.norm(vertices-curr,axis=1)
        norms[i]=1
        norms=np.exp(-(norms**2))
        norms[i]=0
        norms[i]=-np.sum(norms)
        laplacian.append(-norms)
    return np.array(laplacian)
class graph:
    def __init__(self,weights,vertices):
        self.vertices=np.array(vertices)
        self.laplacian=generate_laplacian(vertices)
        self.M=np.diag(weights)
        self.labels=[]
    def cluster(self,k):
        eig_val,eig_vec,=eig(self.laplacian,self.M)
        idx = eig_val.argsort()[:k]
        eig_val=eig_val[idx]
        eig_vec=eig_vec[:,idx]
        eig_vec=normalize(eig_vec,axis=1)
        kmeans = KMeans(n_clusters=k, random_state=0).fit(eig_vec)
        self.labels=kmeans.labels_
def visualize_gaussian():
    mean = (1, 2)
    cov = [[1, 0], [0, 1]]
    x = np.random.multivariate_normal(mean, cov, 50)
    ymean = (6, 5)
    ycov = [[1, 0], [0, 1]]
    y = np.random.multivariate_normal(ymean, ycov, 50)
    gaussian_double=np.concatenate((x,y),axis=0)
    np.random.shuffle(gaussian_double)
    g=graph(50*np.ones(100),gaussian_double)
    g.cluster(2)
    label_zero=gaussian_double[g.labels==0]
    label_one=gaussian_double[g.labels==1]
    axeszero=np.split(np.array(label_zero),2,axis=1)
    zeroxaxis=np.array(axeszero[0]).flatten()
    zeroyaxis=np.array(axeszero[1]).flatten()
    axesone=np.split(np.array(label_one),2,axis=1)
    onexaxis=np.array(axesone[0]).flatten()
    oneyaxis=np.array(axesone[1]).flatten()
    c0 = pl.scatter(zeroxaxis,zeroyaxis,c='r')
    c1 = pl.scatter(onexaxis,oneyaxis,c='g')
    pl.legend([c0, c1], ['gaussian (1,2)', 'gaussian (5,4)'])
    pl.title('gaussian distributions with spectral clustering')
    pl.show()
def visualize_circles():
    radii=[1,6,11]
    circles=[]
    for r in radii:
        circles+=[(math.cos(2*math.pi/50*x)*r,math.sin(2*math.pi/50*x)*r) for x in range(50)]
        circles+=[(math.cos(2*math.pi/50*x)*(r+1),math.sin(2*math.pi/50*x)*(r+1)) for x in range(50)]
    circles=np.array(circles)
    np.random.shuffle(circles)
    g=graph(5*np.ones(len(radii)*100),circles)
    g.cluster(3)
    label_zero=circles[g.labels==0]
    label_one=circles[g.labels==1]
    label_two=circles[g.labels==2]
    axeszero=np.split(np.array(label_zero),2,axis=1)
    zeroxaxis=np.array(axeszero[0]).flatten()
    zeroyaxis=np.array(axeszero[1]).flatten()
    axesone=np.split(np.array(label_one),2,axis=1)
    onexaxis=np.array(axesone[0]).flatten()
    oneyaxis=np.array(axesone[1]).flatten()
    axes2=np.split(np.array(label_two),2,axis=1)
    twoxaxis=np.array(axes2[0]).flatten()
    twoyaxis=np.array(axes2[1]).flatten()
    c0 = pl.scatter(zeroxaxis,zeroyaxis,c='r')
    c1 = pl.scatter(onexaxis,oneyaxis,c='g')
    c2 = pl.scatter(twoxaxis,twoyaxis,c='b')
    pl.legend([c0, c1,c2], ['circle 1', 'circle 2','circle 3'])
    pl.title('circles with spectral clustering')
    pl.show()
def visualize_lines():
    lines=[]
    slopes=[1,2]
    for s in slopes:
        lines+=[(x,s*x) for x in range(50)]
    lines=np.array(lines)
    g=graph(50*np.ones(len(slopes)*50),lines)
    g.cluster(2)
    label_zero=lines[g.labels==0]
    label_one=lines[g.labels==1]
    axeszero=np.split(np.array(label_zero),2,axis=1)
    zeroxaxis=np.array(axeszero[0]).flatten()
    zeroyaxis=np.array(axeszero[1]).flatten()
    axesone=np.split(np.array(label_one),2,axis=1)
    onexaxis=np.array(axesone[0]).flatten()
    oneyaxis=np.array(axesone[1]).flatten()
    c0 = pl.scatter(zeroxaxis,zeroyaxis,c='r')
    c1 = pl.scatter(onexaxis,oneyaxis,c='g')
    pl.legend([c0, c1], ['line 1', 'line 2'])
    pl.title('lines with spectral clustering')
    pl.show()
# visualize_circles()
# visualize_lines()
visualize_gaussian()
