#!/usr/bin/env python3

import awkward as ak
import numpy as np
import os
import shutil
import subprocess

dasgoclient = '/cvmfs/cms.cern.ch/common/dasgoclient'
redirectors = {
    'ucsd': 'root://redirector.t2.ucsd.edu:1095/',
    'ucsd_xcache': 'root://xcache-redirector.t2.ucsd.edu:2042/',
    'fnal': 'root://cmsxrootd.fnal.gov/',
    'fnal_eos': 'root://cmseos.fnal.gov/',
    'global': 'root://cms-xrd-global.cern.ch/',
    'europe': 'root://xrootd-cms.infn.it/'
}

def das_wrapper(DASname, qualifier='dataset', query='file', mask=''):
    sampleName = DASname.rstrip('/')
    dbsOut = subprocess.check_output(
        [dasgoclient, f"-query={query} {qualifier}={sampleName} {mask}"],
        stderr=subprocess.PIPE,
        text=True,
    )
    dbsOut = dbsOut.split('\n')
    if len(dbsOut[-1]) < 1: dbsOut = dbsOut[:-1]
    return dbsOut

def xrdfsls(path, redirector):
    cmd = f"xrdfs {redirector} ls {path}"
    out = os.popen(cmd).readlines()
    out = [ l.replace('\n','') for l in out ]
    return out

def get_four_vec_fromPtEtaPhiM(cand, pt, eta, phi, M, copy=True):
    '''
    Get a LorentzVector from a NanoAOD candidate with custom pt, eta, phi and mass
    All other properties are copied over from the original candidate
    '''
    from coffea.nanoevents.methods import vector
    ak.behavior.update(vector.behavior)

    vec4 = ak.zip(
        {
            "pt": pt,
            "eta": eta,
            "phi": phi,
            "mass": M,
        },
        with_name="PtEtaPhiMLorentzVector",
    )
    if copy:
        vec4.__dict__.update(cand.__dict__)
    return vec4

def delta_phi(first, second):
    return np.arccos(np.cos(first.phi - second.phi))

def delta_r2(first, second):
    return (first.eta - second.eta) ** 2 + delta_phi(first, second) ** 2

def delta_r(first, second):
    return np.sqrt(delta_r2(first, second))

def choose(first, n=2):
    assert n<5, "Can't handle combinations of more than 4 objects"
    tmp = ak.combinations(first, n, fields=["first", "second", "third", "fourth"][:n])
    combs = tmp['0']
    for i in range(1,n):
        combs = combs.__add__(tmp[str(i)])
    for i in range(n):
        combs[str(i)] = tmp[str(i)]
    return combs

def cross(first, second):
    tmp = ak.cartesian([first, second])
    combs = (tmp['0'] + tmp['1'])
    combs['0'] = tmp['0']
    combs['1'] = tmp['1']
    return combs

def match(first, second, deltaRCut=0.4):
    drCut2 = deltaRCut**2
    combs = ak.cartesian([first, second], nested=True)
    return ak.any((delta_r2(combs['0'], combs['1'])<drCut2), axis=2)

def finalizePlotDir( path ):
    path = os.path.expandvars(path)
    if not os.path.isdir(path):
        os.makedirs(path)
    shutil.copy( os.path.expandvars( 'php/index.php' ), path )
