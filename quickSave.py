#Quick save/load library - These functions allow for convenient and fast saving/loading
#of python variables.
#Original author: Joey Wilson

import cPickle

def quickLoad(fileName):
    """Load a variable from pickle file in one line."""
    fid = open(fileName,'r')
    Y = cPickle.load(fid)
    fid.close()
    return Y

def quickSave(fileName,var):
    """Save a variable to a pickle file in one line."""
    fid = open(fileName,'w')
    cPickle.dump(var,fid,2)
    fid.close()
