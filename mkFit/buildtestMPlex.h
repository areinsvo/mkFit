#ifndef _buildtest_mplex_
#define _buildtest_mplex_

#include "Event.h"
#include "Track.h"

/*
void   generateTracks(std::vector<Track>& simtracks, int Ntracks);
void   make_validation_tree(const char         *fname,
                            std::vector<Track> &simtracks,
                            std::vector<Track> &rectracks);
*/

double runBuildingTest(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);
double runBuildingTestBestHit(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);

void buildTestParallel(std::vector<Track>& evt_seeds,std::vector<Track>& evt_track_candidates,
		       std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx);

void processCandidates(std::pair<Track, TrackState>& cand,std::vector<std::pair<Track, TrackState> >& tmp_candidates,
		       unsigned int ilay,std::vector<std::vector<Hit> >& evt_lay_hits,std::vector<std::vector<BinInfo> >& evt_lay_phi_hit_idx);

double runBuildingTestPlex(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);
double runBuildingTestPlexOld(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);
double runBuildingTestPlexBestHit(std::vector<Track>& simtracks/*, std::vector<Track>& rectracks*/);

#endif