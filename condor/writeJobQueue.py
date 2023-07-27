import os,sys

### Data
#DataB: 291987 -> 1.2 with pace=250000
#DataC: 14797731 -> 59.2 with pace=250000
#DataD: 7707588 -> 30.8 with pace=250000
#DataE: 11616645 -> 46.5 with pace=250000
#DataF: 43122993 -> 172.5 with pace=250000
#DataG: 3669414 -> 14.7 with pace=250000

nes = dict()
nes['DataB'] = 291987
nes['DataC'] = 14797731
nes['DataD'] = 7707588
nes['DataE'] = 11616645
nes['DataF'] = 43122993
nes['DataG'] = 3669414

pace = 250000
fout = open("queue.txt","w")
for d in nes.keys():
    tj = 0
    while tj*pace < nes[d]:
        fout.write("$ENV(SCOUTINGINPUTDIR) $ENV(SCOUTINGOUTPUTDIR) --data --partialUnblinding --inSample %s --splitIndex %d\n"%(d,tj)) 
        tj = tj+1
