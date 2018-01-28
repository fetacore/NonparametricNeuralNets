# -*- coding: utf-8 -*-
"""

@author: rook
"""

import numpy as np
import math

#hlayers= a list of numbers of hidden units for arbitrarily many hidden layers
#active_fnct= a list of types of activation functions (one more argument than hlayers)


#Basic auxiliary functions
def comb_arrays(x):
    #x is the list of arrays (for multivariate environment)
    z = np.zeros((len(x[0][:]),len(x)))
    for i in range(0,len(x)):
        for j in range(0,len(x[0])):
            z[j,i] = x[i][j]
    return(z)
        
def w_list_creator(x,hlayers=[2]):
    weights = []
    for i in range(0,(len(hlayers)+1)):
        if i == 0:
            weights.append(np.array(np.random.rand(len(x[0,:]),hlayers[0])))
        elif i == (len(hlayers)):
            weights.append(np.array(np.random.rand(hlayers[len(hlayers)-1],1)))
        else:
            weights.append(np.array(np.random.rand(hlayers[i-1],hlayers[i])))
    return(weights)

def bias_list_creator(hlayers=[2]):
    biases = []
    for i in range(0,(len(hlayers)+1)):
        if i == 0:
            biases.append(np.array(np.random.rand(1,hlayers[0])))
        elif i == (len(hlayers)):
            biases.append(np.array(np.random.rand(1,1)))
        else:
            biases.append(np.array(np.random.rand(1,hlayers[i])))
    return(biases)

def activation(x_i, arg):
    if arg == 'tanh':
        e = math.e
        z = np.zeros((len(x_i[:,0]),len(x_i[0,:])))
        for i in range(0,len(x_i[0,:])):
            for j in range(0,len(x_i[:,0])):
                z[j,i] = (e**(x_i[j,i])-e**(-x_i[j,i]))/(e**(x_i[j,i])+e**(-x_i[j,i]))
    elif arg == 'Sigmoid':
        e = math.e
        z = np.zeros((len(x_i[:,0]),len(x_i[0,:])))
        for i in range(0,len(x_i[0,:])):
            for j in range(0,len(x_i[:,0])):
                z[j,i] = 1/(1+e**(-x_i[j,i]))
    elif arg == 'Softsign':
        z = np.zeros((len(x_i[:,0]),len(x_i[0,:])))
        for i in range(0,len(x_i[0,:])):
            for j in range(0,len(x_i[:,0])):
                z[j,i] = (x_i[j,i])/(1+abs(x_i[j,i]))
    elif arg == 'Softplus':
        e = math.e
        z = np.zeros((len(x_i[:,0]),len(x_i[0,:])))
        for i in range(0,len(x_i[0,:])):
            for j in range(0,len(x_i[:,0])):
                z[j,i] = np.log(1+e**(x_i[j,i]))
    elif arg == 'Arctan':
        z = np.arctan(x_i)
    else:
        z = x_i
    return(z)
    
def gradient(x_i, arg):
    if arg == 'tanh':
        e = math.e
        z = 0
        g = np.zeros((len(x_i[:,0]),len(x_i[0,:])))
        for i in range(0,len(x_i[0,:])):
            for j in range(0,len(x_i[:,0])):
                z = (e**(x_i[j,i])-e**(-x_i[j,i]))/(e**(x_i[j,i])+e**(-x_i[j,i]))
                g[j,i] = 1- z**2
    elif arg == 'Sigmoid':
        e = math.e
        z = 0
        g = np.zeros((len(x_i[:,0]),len(x_i[0,:])))
        for i in range(0,len(x_i[0,:])):
            for j in range(0,len(x_i[:,0])):
                z = 1/(1+e**(-x_i[j,i]))
                g[j,i] = z * (1- z)
    elif arg == 'Softsign':
        g = np.zeros((len(x_i[:,0]),len(x_i[0,:])))
        for i in range(0,len(x_i[0,:])):
            for j in range(0,len(x_i[:,0])):
                g[j,i] = 1/((1+abs(x_i[j,i]))**2)
    elif arg == 'Softplus':
        e = math.e
        g = np.zeros((len(x_i[:,0]),len(x_i[0,:])))
        for i in range(0,len(x_i[0,:])):
            for j in range(0,len(x_i[:,0])):
                g[j,i] = 1/(1+e**(-x_i[j,i]))
    elif arg == 'Arctan':
        g = np.zeros((len(x_i[:,0]),len(x_i[0,:])))
        for i in range(0,len(x_i[0,:])):
            for j in range(0,len(x_i[:,0])):
                g[j,i] = 1/((x_i[j,i])**2 + 1)
    else:
        g = np.ones(x_i.shape)
    return(g)
    
def normalize(x):
    for i in range(0,len(x[0,:])):
        for j in range(0,len(x[:,0])):
            x[j,i] = (x[j,i] - np.mean(x[:,i]))/np.std(x[:,i])
    return(x)


#The real deal
def neurons(x_i,w,b,hlayers=[2],active_fnct=['tanh','linear']):
    #feeds at every node (#hlayers+output layer)
    z = []
    #feeds at every node after activation (#hlayers+output layer)
    active_z = []
    for i in range(0,(len(hlayers)+1)):
        if i == 0:
            z.append(np.dot(x_i,w[i])+b[i])
            active_z.append(activation(z[i],active_fnct[i]))
        elif not(i == (len(hlayers))):
            z.append(np.dot(active_z[i-1],w[i])+b[i])
            active_z.append(activation(z[i],active_fnct[i]))
        else:
            z.append(np.dot(active_z[i-1],w[i])+b[i])
            active_z.append(activation(z[i],active_fnct[len(active_fnct)-1]))
    return(list([z,active_z]))
    
def delta_rule(target,z,w,hlayers=[2],active_fnct=['tanh','linear']):
    #target is the scalar corresponding regressant
    i = len(hlayers)
    deltas = []
    j = 0
    while not(i<0):
        if i == len(hlayers):
            deltas.append(gradient(z[0][i],active_fnct[i])*(z[1][i]-target))
        elif i < 0:
            break
        else:
            o = gradient(z[0][i],active_fnct[i])
            o = o.T
            e = np.dot(w[i+1],deltas[j-1])
            deltas.append(e*o)
        i = i - 1
        j = j + 1
        
    return(deltas)

def update(x,z,w,b,deltas,eta=0.005,hlayers=[2]):
    new_w = []
    new_b = []
    j = len(hlayers)
    for i in range(0,len(hlayers)+1):
        if i == 0:
            gw = np.dot(deltas[j],x)
        else:
            gw = np.dot(deltas[j],z[1][i-1])
        new_w.append(w[i]-eta*gw.T)
        gb = deltas[j]
        new_b.append(b[i]-eta*gb.T)
        j = j - 1
    return(list([new_w,new_b]))

def pred_loss(target,z):
    loss = 0.5 * (np.sum(z-target))**2
    return(loss)

def extraction(pred,loss):
    loss = np.min(loss)
    weights = pred[loss][1]
    biases = pred[loss][2]
    return([weights,biases])

def predict(new_x,w,b,hlayers,active_fnct):
    predictor = np.zeros((len(new_x[:,0]),1))
    for i in range(0,len(new_x[:,0])):
        x_i = new_x[i,:]
        x_i.shape = (1,len(x_i))
        z=neurons(x_i,w,b,hlayers,active_fnct)
        predictor[i] = z[1][:]
    return(predictor)

#The Model

##Serial Correlation Case
#n = 100
#x_1 = np.random.normal(5,3,size=(n,1))
#x_2 = 5 * x_1 + np.random.normal(2,1,size=(n,1))
#x = comb_arrays([x_1,x_2])
#target = x_1**2 + x_2 + np.random.normal(0,2,size=(n,1))
#eta = 0.005
#iterations=400
#hlayers = [10,50,500,3]
#active_fnct = ['tanh','tanh','Sigmoid','tanh','linear']

##Serial Correlation and extra variables
#n = 100
#x_1 = np.random.normal(5,3,size=(n,1))
#x_2 = 5 * x_1 + np.random.normal(2,1,size=(n,1))
#x_3 = np.random.normal(0,3,size=(n,1))
#x = comb_arrays([x_1,x_2,x_3])
#target = x_1**2 + x_2 + np.random.normal(0,2,size=(n,1))
#eta = 0.005
#iterations=100
#hlayers = [10,50,500,3]
#active_fnct = ['tanh','tanh','Sigmoid','tanh','linear']

##Ommitted Variables and Irrelevant
#n = 100
#x_1 = np.random.normal(5,3,size=(n,1))
#x_2 = 5 * x_1 + np.random.normal(2,1,size=(n,1))
#x_3 = np.random.normal(0,3,size=(n,1))
#x = comb_arrays([x_1,x_3])
#target = x_1**2 + x_2 + np.random.normal(0,2,size=(n,1))
#eta = 0.005
#iterations=100
#hlayers = [10,50,500,3]
#active_fnct = ['tanh','tanh','Sigmoid','tanh','linear']

#Ommitted Variables
#n = 100
#x_1 = np.random.normal(5,3,size=(n,1))
#x_2 = 5 * x_1 + np.random.normal(2,1,size=(n,1))
#x_3 = np.random.normal(0,3,size=(n,1))
#x = x_2
#target = x_1**2 + x_2 + np.random.normal(0,2,size=(n,1))
#eta = 0.005
#iterations=100
#hlayers = [10,50,500,3]
#active_fnct = ['tanh','tanh','Sigmoid','tanh','linear']


#n = 100
#x_1 = np.random.normal(5,3,size=(n,1))
#x_2 = np.random.normal(1,3,size=(n,1))
#true_model = x_1**2 + np.sign(x_2)
#error = np.random.normal(0,1,size=(n,1))
#target = true_model + error
#x = comb_arrays([x_1,x_2])
#n = 1000
#x = np.random.normal(5,3,size=(n,1))
##x_2 = np.random.normal(1,3,size=(n,1))
#true_model = x**2
#error = np.random.normal(0,5,size=(n,1))
#target = true_model + error
#threshold = 1e-24
#eta = 0.005
#hlayers = [100,50,50,30,20,5,30,2]
#active_fnct = ['tanh','Arctan','Softsign','Sigmoid','tanh','tanh','tanh','tanh','linear']
#hlayers = [50,20]
#active_fnct = ['tanh','tanh','linear']
#[zreg,w,b,loss,iterations] = neural_sifu(target,x,threshold,eta,hlayers,active_fnct)
#pred = perceptron(x[len(x)-1,:],w,b,hlayers,active_fnct)
#i = np.arange(0,iterations,1)
#plt.figure(1)
#plt.plot(i,loss)
#plt.xlabel('Iterations')
#plt.ylabel('Loss')
#plt.title('Loss from data')
#plt.grid(True)
#plt.show('hold')
