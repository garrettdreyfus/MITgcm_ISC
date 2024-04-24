import os, pickle, glob
from scipy.io import loadmat
import numpy as np
from xmitgcm import open_mdsdataset
from pathlib import Path
import re



def matVarsFile(fname):
    one_up = os.path.dirname(fname)
    return one_up+"/input/metaparameters.mat"

def grabMatVars(fname,val_tup):
    variables = loadmat(matVarsFile(fname),variable_names=val_tup,struct_as_record=False)
    return variables

    
def outPath(resultspath):
    nameparts = resultspath.split("/")
    shortname = nameparts[-3] + "|" + nameparts[-2]
    fpath = "/".join((nameparts[:-5]+["pics"] + [shortname]))
    return shortname, fpath

def getIterNums(fpath):
    iters = []
    saltiters = []
    for fname in glob.glob(fpath+"/*.data"):
        if 'THETA' in fname and "inst" not in fname:
            n = int(fname.split(".")[1])
            iters.append(n)
        if 'SALT' in fname and "inst" not in fname:
            n = int(fname.split(".")[1])
            saltiters.append(n)
    return np.unique(np.asarray(np.intersect1d(iters,saltiters)))[:-1]


def grabDeltaT(fname):
    p = Path(fname).parents[0]
    p = str(p)+"/input/data"
    with open(p,"rb") as file:
        for line in file:
            s = str(line)
            if "deltaT" in s:
                result = re.search('=(.*),', s)
                result = result.group(1)
                if result == "0.25000+02":
                    result = "0.25000e+02"
                return float(result)
       

