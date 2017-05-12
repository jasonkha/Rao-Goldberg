import numpy as np
import matplotlib.pyplot as plt
import math
import random
#function sin(x)+cos(5x)+sin(3x)
#we picked a random starting point and used the three methods
#repeated this x*1000 times and saw which method ended up at the lowest spot each times
#kept count for x in range(50)
def func(x):
    return math.sin(x)+math.cos(5*x)+math.sin(3*x)
def accelerated(x,iterations):
    epsilon=.8
    mu=.8
    deltax_prev=0
    for j in range(iterations):
        deltax = mu* deltax_prev - epsilon * (-5*math.sin(5*x)+3*math.cos(3*x)+math.cos(x))
        x_nesterov = deltax + mu* (deltax -deltax_prev) #nesterov update
        x += x_nesterov
        deltax_prev=deltax
        epsilon=.8/(j+1)
    return func(x)

def momentum(x,iterations):
    epsilon=.8
    mu=.8
    deltax_prev=0
    for j in range(iterations):
        #mu*Deltax_prev gives us "momentum"
        deltax = mu* deltax_prev - epsilon * (-5*math.sin(5*x)+3*math.cos(3*x)+math.cos(x))
        x += deltax
        deltax_prev=deltax
        epsilon=.8/(j+1)
    return func(x)
def normal(x,iterations):
    epsilon=.8
    for j in range(iterations):
        deltax =  - epsilon * (-5*math.sin(5*x)+3*math.cos(3*x)+math.cos(x))
        x += deltax
        epsilon=.8/(j+1)
    return func(x)
#plot/iterations
xaxis=[]
a=[]
m=[]
n=[]
for j in range(50):
    counts=[0,0,0]
    for i in range(1000*j):
        x=random.uniform(-50, 50)
        win=np.argmin([accelerated(x,100),momentum(x,100),normal(x,100)])
        counts[win]+=1
    xaxis.append(j*1000)
    a.append(counts[0])
    m.append(counts[1])
    n.append(counts[2])
plt.figure()
plt.plot(xaxis,a, label="accelerated")
plt.plot(xaxis,m, label="momentum")
plt.plot(xaxis,n, label="normal")
plt.legend()
plt.xlabel('number of rounds of random starting points')
plt.ylabel('number of rounds ending up with lowest value')
plt.title('descent method comparisons')
plt.show()
