from analysis import *
import numpy as np
from chooseSurfaces import *
import copy
from scipy.optimize import fmin

# this file is to fit weighted indices also with alternate residuals

def altCost(p,x,y,c,npower):
    # this is to not use nonlinear least squares but just find a 
    return np.sum((abs((fitFunction(x,p,correction=c)-y)/y))**npower)

def fitWeightedData(fileName,indices_to_weight,p0=np.array([0.01,0.02,0.01,1.0]),correction=1,weight=50,altnpower=4,lsq=False):
    indices, unrel_es, rel_es, sizes = loadFileIntoLists(fileName)
    new_indices, new_energies = pruneIndicesAndEnergies(indices,rel_es)
    weighted_indices = copy.copy(new_indices)
    weighted_energies = copy.copy(new_energies)     
    for miller in indices_to_weight:
        ind = new_indices.index(miller)
        ener = new_energies[ind]
        weighted_indices += weight*[miller]
        weighted_energies += weight*([ener])
    if lsq:
	nterms = int(len(p0)-correction)
	print nterms
        output = fitSurfaceEnergies(np.array(weighted_indices),np.array(weighted_energies),n=nterms,p0=p0,correction=correction)
        bfparams = output[0]
        cost_per_index = output[2] 
    else:
        bfparams = fmin(altCost,p0,args=(np.array(weighted_indices),np.array(weighted_energies),correction,altnpower))
        cost_per_index = altCost(bfparams,np.array(weighted_indices),np.array(weighted_energies),correction,altnpower)/len(weighted_indices) 
    plotindices,plotenergies = expandList(new_indices,new_energies)
    plotSubSet(plotindices,plotenergies,[1,-1,0],[1,1,0],bfparams,correction=correction)
    plotSubSet(plotindices,plotenergies,[1,-1,2],[1,1,0],bfparams,correction=correction)
    plotSubSet(plotindices,plotenergies,[1,1,-1],[0,1,1],bfparams,correction=correction)
    
    return bfparams, cost_per_index 
    
