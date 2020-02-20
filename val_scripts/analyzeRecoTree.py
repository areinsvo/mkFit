import os, sys, ROOT, math
#Run with /usr/bin/python2.7 analyzeRecoTree.py

sw = ROOT.TStopwatch()
sw.Start()


chain_in = ROOT.TChain("recotree")
chain_in.Add("valtree_1000.root")
out_file = ROOT.TFile("Duplicates_post.root", "recreate")

starting_event =           0000
process_n_events =      1000
if (starting_event + process_n_events) > chain_in.GetEntries():
    last_event = chain_in.GetEntries()
else:
    last_event = starting_event + process_n_events

print "Going to process from ",starting_event," to ",last_event

h_nHitsMatched_duplicate = ROOT.TH1F("nHitsMatched_duplicate","",30,0,30)
h_fracHitsMatched_duplicate = ROOT.TH1F("fracHitsMatched_duplicate","",100,0,1)
h_nHitsMatched = ROOT.TH1F("nHitsMatched","",30,0,30)
h_fracHitsMatched = ROOT.TH1F("fracHitsMatched","",50,0,1)
h_nHits_duplicate = ROOT.TH1F("nHits_duplicate","",30,0,30)
h_nHits = ROOT.TH1F("nHits","",30,0,30)
h_hitlyrs = ROOT.TH1F("hitlyrs","",100,0,100)

h_dPt_duplicates = ROOT.TH1F("dPt_duplicates","",100,0,10)
h_frac_dPt_duplicates = ROOT.TH1F("frac_dPt_duplicates","",100,0,1)

h_dPt = ROOT.TH1F("dPt","",80,0,80)
h_frac_dPt = ROOT.TH1F("frac_dPt","",100,0,1)

h_dPhi = ROOT.TH1F("dPhi","",100,0,0.5)
h_dPhi_duplicates = ROOT.TH1F("dPhi_duplicates","",100,0,0.5)

h_eta_duplicate = ROOT.TH1F("eta_duplicate","",250,0,2.5)

h_dEta = ROOT.TH1F("dEta","",200,0,2)
h_dEta_duplicates = ROOT.TH1F("dEta_duplicates","",200,0,2)

h_nHitsShared_duplicates = ROOT.TH1F("nHitsShared_duplicates","",30,0,30)
h_fracHitsShared_duplicates = ROOT.TH1F("fracHitsShared_duplicates","",100,0,1)

h_nHitsShared = ROOT.TH1F("nHitsShared","",30,0,30)
h_fracHitsShared = ROOT.TH1F("fracHitsShared","",100,0,1)

total_Duplicates = 0;
total_Duplicates2 = 0;
total_Duplicates_Removed = 0;
total_nonDuplicates = 0;
total_nonDuplicates_Removed = 0;
total_RecoTracks = 0;

maxdEta = 0.2 #default 0.2
maxdPhi = 0.2 #default 0.1
maxdPt = 1.0 #default 0.05
minFracSharedHits = 0.75

#Loop over events
for j_entry in range(starting_event,last_event):
    i_entry = chain_in.LoadTree(j_entry)

    if i_entry < 0:
        break

    nb = chain_in.GetEntry(j_entry)
    if nb <= 0:
        continue
    
    if j_entry % 10 ==0:
        print "Processing entry " + str(j_entry) + " of " + str(last_event)
    for i in range(0,chain_in.mcID.size()):
        total_RecoTracks += 1
        removed=False
        if chain_in.isDuplicate[i]:
            h_nHitsMatched_duplicate.Fill(chain_in.nHitsMatched[i])
            h_fracHitsMatched_duplicate.Fill(chain_in.fracHitsMatched[i])
            h_nHits_duplicate.Fill(chain_in.nHits[i])
            h_eta_duplicate.Fill(abs(chain_in.eta[i]))
            total_Duplicates += 1
        else:
            h_nHitsMatched.Fill(chain_in.nHitsMatched[i])
            h_fracHitsMatched.Fill(chain_in.fracHitsMatched[i])
            h_nHits.Fill(chain_in.nHits[i])
            total_nonDuplicates += 1

        mc1 = chain_in.mcID[i]
        pt1 = chain_in.pt[i]
        phi1 = chain_in.phi[i]
        minPtDiff = 10000
        minPhiDiff = 10
        minEtaDiff = 10
        eta1 = chain_in.eta[i]

        for j in range(i+1,chain_in.mcID.size()):
            mc2 = chain_in.mcID[j]
            pt2 = chain_in.pt[j]
            phi2 = chain_in.phi[j]
            eta2 = chain_in.eta[j]
            dphi = abs(math.acos(math.cos(phi2-phi1)))
            deta = abs(eta2-eta1)
            duplicatePreselection = (dphi < maxdPhi and abs(pt2-pt1)/max(pt1,pt2) < maxdPt and deta < maxdEta) #and abs(chain_in.eta[i]) >1.4)
#duplicatePreselection = (dphi < maxdPhi and deta < maxdEta) #and abs(chain_in.eta[i]) >1.4)

            if mc1 == mc2:
                if not chain_in.isDuplicate[i] or not chain_in.isDuplicate[j]:
                    print "ERROR!"
                
                total_Duplicates2 += 1
                h_dPt_duplicates.Fill(abs(pt2-pt1))
                h_frac_dPt_duplicates.Fill(abs(pt2-pt1)/pt1)
                h_dPhi_duplicates.Fill(dphi)
                h_dEta_duplicates.Fill(deta)
                hit_list = []
                numHits = 0.0
                numHitsShared = 0.0
                #print("Here!")
                for k in range(chain_in.hitidxs[i].size()):
                    hit_list.append(chain_in.hitidxs[i][k])
                    numHits += 1
                for m in range(chain_in.hitidxs[j].size()):
                    if chain_in.hitidxs[j][m] in hit_list:
                        numHitsShared+= 1
                h_nHitsShared_duplicates.Fill(numHitsShared)
                h_fracHitsShared_duplicates.Fill(numHitsShared/numHits)
                
                fracHitsShared = numHitsShared/numHits;
                
                if duplicatePreselection and fracHitsShared > minFracSharedHits:
                    total_Duplicates_Removed += 1
            elif not duplicatePreselection:
                continue;
            else:
                minPtDiff = min(minPtDiff,abs(pt2-pt1))
                minPhiDiff = min(minPhiDiff,dphi)
                minEtaDiff = min(minEtaDiff,deta)
                numHits = 0.0
                numHitsShared = 0.0
                hit_list = []

                for k in range(chain_in.hitidxs[i].size()):
                    hit_list.append(chain_in.hitidxs[i][k])
                    numHits+= 1
                for m in range(chain_in.hitidxs[j].size()):
                    if chain_in.hitidxs[j][m] in hit_list:
                        numHitsShared+= 1
                h_nHitsShared.Fill(numHitsShared)
                h_fracHitsShared.Fill(numHitsShared/numHits)

                fracHitsShared = numHitsShared/numHits
                if fracHitsShared > minFracSharedHits:
                    removed = True
                    
        if removed and not chain_in.isDuplicate[i]:
            total_nonDuplicates_Removed += 1


        h_dPt.Fill(minPtDiff)
        h_frac_dPt.Fill(minPtDiff/pt1)
        h_dPhi.Fill(minPhiDiff)
        h_dEta.Fill(minEtaDiff)
            #else:
#   print "Checking if it is close"
print "Parameters: maxdEta = "+ str(maxdEta) + ", maxdPhi = " + str(maxdPhi)+ ", maxdPt = " + str(maxdPt) + ", minFracSharedHits = " + str(minFracSharedHits)
print "Total number of duplicates is ", total_Duplicates, " or ", total_Duplicates2
print "Total number of duplicates removed is ", total_Duplicates_Removed
print "Total number of non-duplicates is ", total_nonDuplicates
print "Total number of non-duplicates removed is ", total_nonDuplicates_Removed
print "Total number of reconstructed tracks is ", total_RecoTracks

print str(maxdEta) + " " + str(maxdPhi)+ " " + str(maxdPt) + " " + str(minFracSharedHits) + " " + str(total_Duplicates_Removed)  + " " + str(total_nonDuplicates_Removed) + " " + str(total_Duplicates) + " " + str(total_RecoTracks)

#for j in range(0,chain_in.hitlyrs[i].size()):
#h_hitlyrs.Fill(chain_in.hitlyrs[i][j])



        #if chain_in.duplmask_build < 0:
        #continue
    

out_file.Write()
out_file.Close()

sw.Stop()

#print automatically appends a new line character
print "Real Time = " + str(sw.RealTime() / 60.0 ) + " minutes."
print "CPU Time = " + str(sw.CpuTime() / 60.0 ) + " minutes."
