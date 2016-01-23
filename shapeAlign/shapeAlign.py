from __future__ import division
import pandas as pd
import itertools
import scipy
import mahotas
import sklearn
import seaborn as sns
import glob
import os.path as op
import os

def getSubMatrix(A,B,s):
    """Return linearized submatrices of the linearized matrices A and B.

    Parameters
    ----------
    A : pd.core.frame.DataFrame
        linearized matrix of dimensions (1xn)
    B : pd.core.frame.DataFrame
        linearized matrix of dimensions (1xm)
    s : int
        shift of A relative to B

    Returns
    -------
    Asub : np.ndarray
        linearized submatrix of A
    Bsub : np.ndarray
        linearized submatrix of B
    BsubR : np.ndarray
        linearized, reversed submatrix of B
    """

    # Convert to arrays and replace NaNs with zero
    arrA = np.nan_to_num(np.array(A))
    arrB = np.nan_to_num(np.array(B))

    # Reshape
    matA = np.reshape(arrA,(4,200))
    matB = np.reshape(arrB,(4,200))
    matB_r = np.fliplr(matB)

    if (s<0):
        Asub = matA[0:,0:s]
        Bsub = matB[0:,abs(s):]
        BsubR = matB_r[0:,abs(s):]
    elif (s>0):
        Asub = matA[0:,0:-s]
        Bsub = matB[0:,0:200-s]
        BsubR = matB_r[0:,0:200-s]
    else:
        Asub = matA[0:,0:]
        Bsub = matB[0:,0:]
        BsubR = matB_r[0:,0:]

    return np.ravel(Asub),np.ravel(Bsub),np.ravel(BsubR)

def getOptShift(A,B,smin,smax):
    """Return the optimal shift for two 1D representations of 2D arrays

    Parameters
    ----------
    A : pd.core.frame.DataFrame
        linearized matrix of dimensions (1xn)
    B : pd.core.frame.DataFrame
        linearized matrix of dimensions (1xm)
    smin : int
        minimum shift
    smax : int
        maximum shift

    Returns
    -------
    optShift : int
        Shift associated with minimum distance
    optDist : double
        Minimum distance
    optRev : bool
        Is the optimal alignment reversed?
    """

    shifts=np.arange(smin,smax+1)
    dist=[]
    rev=[]
    for s in shifts:
        Asub,Bsub, BsubR = getSubMatrix(A,B,s)
        fDist = scipy.spatial.distance.correlation(Asub,Bsub)
        rDist = scipy.spatial.distance.correlation(Asub,BsubR)
        if (fDist >= rDist):
            dist.append(fDist)
            rev.append(False)
        else:
            dist.append(rDist)
            rev.append(True)

    optidx = np.argmin(dist)
    optShift = shifts[optidx]
    optDist = dist[optidx]
    optRev = rev[optidx]
    return optShift, optDist, optRev


def processShapeFile(Fn):
    """
    Process a given shape features file (each row is a site, and each column is a shape parameter).
    Note that this file should contain all of the shape parameters as a vector. It is assumed that
    the window length is 200. This should be fixed in future versions to allow differing window widths.

    Parameters
    ----------
    Fn : str
        Filename

    Returns
    -------
    siteList : list
        List of sites; order corresponds to row/column ordering of matrices that are also returned
    distMat : np.ndarray
        Distance matrix summarizing pairwise distances between sites
    shiftMat : np.ndarray
        Matrix summarizing optimal shifts between sites
    revMat : np.ndarray
        Matrix summarizing whether alignment between two sites required reversal
    """

    shapeDf = pd.read_csv(Fn,header=None,sep='\t')
    siteList=shapeDf[0]
    distMat = np.zeros((len(siteList),len(siteList)))
    shiftMat = np.zeros((len(siteList),len(siteList)))
    revMat = np.zeros((len(siteList),len(siteList)),dtype=bool)

    for i in range(0,len(siteList)):
        print i, '\t', siteList[i]
        for j in range(i,len(siteList)):
            shiftMat[i,j], distMat[i,j], revMat[i,j] = getOptShift(shapeDf.loc[shapeDf[0]==siteList[i]].ix[:,1:],shapeDf.loc[shapeDf[0]==siteList[j]].ix[:,1:],-25,25)
            shiftMat[j,i] = shiftMat[i,j]
            distMat[j,i] = distMat[i,j]
            revMat[j,i] = revMat[i,j]

    return siteList, distMat, shiftMat, revMat

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', help='Name of the directory to be processed', type=str)
    args = parser.parse_args()
    processFolder(args.folder)
