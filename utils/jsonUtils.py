import os,sys,json

def isgoodrun(trun,tlumiBlock,runs):
    srun="%d"%trun
    if srun not in runs.keys():
        return False
    for lumiBlock in runs[srun]:
        if len(lumiBlock) != 2 :
            print('ERROR reading lumiBlock from JSON: run:',run,'lumiBlock:',lumiBlock)
            exit()
        if tlumiBlock >= int(lumiBlock[0]) and tlumiBlock <= int(lumiBlock[1]):
            return True
    return False
