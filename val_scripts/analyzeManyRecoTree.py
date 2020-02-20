import os, sys, ROOT, math
#Run with /usr/bin/python2.7 analyzeRecoTree.py

sw = ROOT.TStopwatch()
sw.Start()


chain_in = ROOT.TChain("recotree")
chain_in.Add("valtree_1000.root")
out_file = ROOT.TFile("Duplicates_many.root", "recreate")

starting_event =        0000
process_n_events =      1000
if (starting_event + process_n_events) > chain_in.GetEntries():
    last_event = chain_in.GetEntries()
else:
    last_event = starting_event + process_n_events

print "Going to process from ",starting_event," to ",last_event

total_Duplicates = 0;
total_Duplicates_Removed = 0;
total_nonDuplicates_Removed = 0;
total_RecoTracks = 0;

#maxdEta = [0.05,0.1,0.2] #default 0.2
#maxdPhi = [0.2]#default 0.1
#maxdPt = [0.25,0.3,0.4] #default 0.05
minFracSharedHits = [0.75] #default 0.75
maxdR = 0.003
maxdEta = [0.1] #default 0.2
maxdPhi = [0.4]#default 0.1
maxdPt = [0.7] #default 0.05
#minFracSharedHits = [0.75] #default 0.75

#maxdEta = [0.05,0.1,0.15,0.2,0.25] #default 0.2
#maxdPhi = [0.4,0.425,0.45,0.5]#default 0.1
#maxdPt = [0.7] #default 0.05
#minFracSharedHits = [0.75] #default 0.75

dup_removed = [[[[0 for t in range(0,len(minFracSharedHits))]
                  for z in range(0,len(maxdPt))]
                 for y in range(0,len(maxdPhi))]
                for x in range(0,len(maxdEta))]

good_non_dup_removed = [[[[0 for t in range(0,len(minFracSharedHits))]
                      for z in range(0,len(maxdPt))]
                    for y in range(0,len(maxdPhi))]
                   for x in range(0,len(maxdEta))]

fakes_removed = [[[[0 for t in range(0,len(minFracSharedHits))]
                     for z in range(0,len(maxdPt))]
                    for y in range(0,len(maxdPhi))]
                   for x in range(0,len(maxdEta))]

numPairsHitCompare = [[[[0 for t in range(0,len(minFracSharedHits))]
                   for z in range(0,len(maxdPt))]
                  for y in range(0,len(maxdPhi))]
                 for x in range(0,len(maxdEta))]

#Loop over events
for j_entry in range(starting_event,last_event):
    i_entry = chain_in.LoadTree(j_entry)

    if i_entry < 0:
        break

    nb = chain_in.GetEntry(j_entry)
    if nb <= 0:
        continue
    
#if j_entry != 558: continue
    if j_entry % 1 ==0:
        print "Processing entry " + str(j_entry) + " of " + str(last_event)
    
    #removed[chain_in.mcID.size()] = false
#removed = []
    removed = [[[[[False for t in range(0,len(minFracSharedHits))]
                 for z in range(0,len(maxdPt))]
                 for y in range(0,len(maxdPhi))]
                    for x in range(0,len(maxdEta))]
                        for i in range(0,chain_in.mcID.size())]

        #removed = [[[ [False for z in range(0,chain_in.mcID.size())]
        #             for y in range(0,len(maxdEta))]
        # for x in range(0,len(maxdPhi))]
#           for i in range(0,len(maxdPt))]

#removed [ track num ] [maxdeta] [maxdphi] [max dpt]

#removed[chain_in.mcID.size()][maxdEta.size()][maxdPhi.size()][ maxdPt.size()] = False

    #loop over built tracks
    for i in range(0,chain_in.mcID.size()):
        if chain_in.mcID[i] < 0 and not (chain_in.mcID[i] == -9 or chain_in.mcID[i] == -1 ): continue
        # --------------------------------
        total_RecoTracks += 1
        # --------------------------------

        #   removed=False
        
        if chain_in.isDuplicate[i] and chain_in.mcID[i] > 0:
            # --------------------------------
            total_Duplicates += 1
            # --------------------------------

        mc1 = chain_in.mcID[i]

        phi1 = chain_in.phi[i]
        eta1 = chain_in.eta[i]
        pt1 = chain_in.pt[i]

        #print "Hello world"
        
        #inner loop over the rest of the built tracks
        for j in range(i+1,chain_in.mcID.size()):
            if chain_in.mcID[j] < 0 and not (chain_in.mcID[j] == -9 or chain_in.mcID[j] == -1 ): continue


            mc2 = chain_in.mcID[j]
            pt2 = chain_in.pt[j]
            phi2 = chain_in.phi[j]
            eta2 = chain_in.eta[j]
            dphi = abs(math.acos(math.cos(phi2-phi1)))
            if dphi > 3.1415926: print "Maybe phi is the problem", str(chain_in.seedID[i])," ", str(chain_in.seedID[j])," ",dphi
            deta = abs(eta2-eta1)
            dFracPt = abs((pt2-pt1)/max(pt1,pt2))
            
            dR = math.sqrt(math.pow(deta,2) + math.pow(dphi,2))


            
            weird = False
            if chain_in.seedID[i] == 153 and chain_in.seedID[j] == -518:
                print phi1, phi2, dphi, maxdPhi[y]
                weird = True


            
            for x in range(0,len(maxdEta)):
                for y in range(0,len(maxdPhi)):
                    for z in range(0,len(maxdPt)):
                        if dR < maxdR:
                            removed[i][x][y][z][t] = True
                            removed[j][x][y][z][t] = True
                            continue
                        
                        if (dphi < maxdPhi[y] and dFracPt < maxdPt[z] and deta < maxdEta[x]) :
                            numPairsHitCompare[x][y][z][t] += 1
                            #print "Hello!" + str(maxdEta[x]) + " " + str(maxdPhi[y])+ " " + str(maxdPt[z])
                            hit_list = []
                            numHits1 = 0.0
                            numHitsShared = 0.0
                            numHits2 = 0.0
                            if weird:
                                print "Hits in track 1:"
                            for k in range(chain_in.hitidxs[i].size()):
                                if chain_in.hitidxs[i][k] >= 0:
                                    #Need both the index and layer to check if two hits match
                                    pair = (chain_in.hitidxs[i][k], chain_in.hitlyrs[i][k])
                                    hit_list.append(pair)
                                    numHits1 += 1
                                    if weird:
                                        print pair
                                            #   if numHits1 != chain_in.nHits[i]: print "Maybe nHits is the problem: ",str(numHits1)," ",str(chain_in.nHits[i])
                            if weird:
                                print "Hits in track 2:"
                            for m in range(chain_in.hitidxs[j].size()):
                                if chain_in.hitidxs[j][m] >= 0:
                                    numHits2 += 1
                                    pair = (chain_in.hitidxs[j][m], chain_in.hitlyrs[j][m])

                                    if weird:
                                        print pair
                                    if pair in hit_list:
                                        numHitsShared+= 1
                            if weird:
                                print pt1, pt2, phi1, phi2, eta1, eta2
                                print "Maybe nHits is the problem: ",str(numHits1)," ",str(chain_in.nHits[i])
                                print str(numHitsShared)," ",str(numHits1)," ",str(numHits2)," ",str(chain_in.nHits[j])
                            for t in range(0,len(minFracSharedHits)):
                                if(numHitsShared/min(numHits1,numHits2) > minFracSharedHits[t]):
                                    removed[i][x][y][z][t] = True
                                    removed[j][x][y][z][t] = True
                                        #if (chain_in.isDuplicate[i] or chain_in.isDuplicate[j]):
#print "Removing: " + str(chain_in.seedID[i]) + " " + str(chain_in.seedID[j])
                                    

    for x in range(0,len(maxdEta)):
        for y in range(0,len(maxdPhi)):
            for z in range(0,len(maxdPt)):
                for t in range(0,len(minFracSharedHits)):
                    for i in range(0,chain_in.mcID.size()):
                        if removed[i][x][y][z][t]:
                            if chain_in.mcID[i] > 0:
                                if chain_in.isDuplicate[i]:
                                    dup_removed[x][y][z][t] += 1;
                                #print "Removed", str(chain_in.seedID[i])
                                else:
                                    good_non_dup_removed[x][y][z][t] += 1;
                            else:
                                fakes_removed [x][y][z][t] += 1;


print "maxdEta maxdPhi maxdPt minFracSharedHits dup_removed good_non_dup_removed  fakes_removed total_Duplicates total_RecoTracks"

for x in range(0,len(maxdEta)):
    for y in range(0,len(maxdPhi)):
        for z in range(0,len(maxdPt)):
            for t in range(0,len(minFracSharedHits)):
                print str(maxdR) + " " + str(maxdEta[x]) + " " + str(maxdPhi[y])+ " " + str(maxdPt[z]) + " " + str(minFracSharedHits[t]) + " " + str(dup_removed[x][y][z][t])  + " " + str(good_non_dup_removed[x][y][z][t]) + " " + str(fakes_removed[x][y][z][t]) + " " + str(total_Duplicates) + " " + str(total_RecoTracks) + " " + str(numPairsHitCompare[x][y][z][t])


#print "Parameters: maxdEta = "+ str(maxdEta) + ", maxdPhi = " + str(maxdPhi)+ ", maxdPt = " + str(maxdPt) + ", minFracSharedHits = " + str(minFracSharedHits)
#print "Total number of duplicates is ", total_Duplicates
#print "Total number of duplicates removed is ", total_Duplicates_Removed
#print "Total number of non-duplicates is ", total_nonDuplicates
#print "Total number of non-duplicates removed is ", total_nonDuplicates_Removed
#print "Total number of reconstructed tracks is ", total_RecoTracks

#print str(maxdEta) + " " + str(maxdPhi)+ " " + str(maxdPt) + " " + str(minFracSharedHits) + " " + str(total_Duplicates_Removed)  + " " + str(total_nonDuplicates_Removed) + " " + str(total_Duplicates) + " " + str(total_RecoTracks)
    

out_file.Write()
out_file.Close()

sw.Stop()

#print automatically appends a new line character
print "Real Time = " + str(sw.RealTime() / 60.0 ) + " minutes."
print "CPU Time = " + str(sw.CpuTime() / 60.0 ) + " minutes."
