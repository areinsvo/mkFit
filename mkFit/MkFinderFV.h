#ifndef MkFinderFV_h
#define MkFinderFV_h

#include "MkBase.h"

#include "Track.h"

class CandCloner;
class CombCandidate;
class LayerOfHits;
class FindingFoos;

//#include "Event.h"

//#include "HitStructures.h"

// NOTES from MkFitter ... where things were getting super messy.
//
// Would like to get rid of msErr, msPar arrays ... single one should do.
// Would also like to get rid of Nhits ... but this is a tougher one.
// Fitting is somewhat depending on this (not the new inter-slurp version).
// No need to have hit array as matriplex. Ideally, don't even want the full
// array if we could copy / clone it from the original candidate.
// Copy as needed or copy full size (compile time constant)? memcpy, std::copy 
// Do need per track size.
//
// Can also have nFoundHits ... then do not need countValid/InvalidHits().
//
// For now will prefix the new variables with 'tf' for track finding.
// Changes will go into all finding functions + import-with-hit-idcs.
//
// Actually ... am tempted to make MkFinder :)

template<int nseeds, int ncands>
class alignas(64) MkFinderFV : public MkBase
{
public:
  static constexpr int NNFV = nseeds * ncands;
  static constexpr int Seeds = nseeds;
  static constexpr int Candidates = ncands;

  static constexpr int MPlexHitIdxMax = 16;

  using MPlexHitIdx = Matriplex::Matriplex<int, MPlexHitIdxMax, 1, NNFV>;
  using MPlexQHoT   = Matriplex::Matriplex<HitOnTrack, 1, 1, NNFV>;

  using MPlexHV     = Matriplex::Matriplex<float, HH,  1, NNFV>;
  using MPlexHS     = Matriplex::MatriplexSym<float, HH,  NNFV>;
  using MPlexQF     = Matriplex::Matriplex<float, 1, 1, NNFV>;
  using MPlexQI     = Matriplex::Matriplex<int,   1, 1, NNFV>;

  struct IdxChi2List
  {
    IdxChi2List() {}
    IdxChi2List(int t, int h, int n, float c) : trkIdx(t), hitIdx(h), nhits(n), chi2(c) {}
    int   trkIdx; // candidate index
    int   hitIdx; // hit index
    int   nhits;  // number of hits (used for sorting)
    float chi2;   // total chi2 (used for sorting)
  };

  //----------------------------------------------------------------------------

  MPlexQF    Chi2;
  MPlexQI    Label;   // seed index in global seed vector (for MC truth match)

  MPlexQI    NHits;
  MPlexQI    NFoundHits;
  HitOnTrack HoTArrs[NNFV][Config::nMaxTrkHits];

  MPlexQI    SeedIdx; // seed index in local thread (for bookkeeping at thread level)
  MPlexQI    CandIdx; // candidate index for the given seed (for bookkeeping of clone engine)

  // Hit indices into LayerOfHits to explore.
  MPlexQI     XHitSize;
  MPlexHitIdx XHitArr;
  MPlexQF     XHitChi2[MPlexHitIdxMax];

  // Hit error / parameters for hit matching, update.
  MPlexHS    msErr;
  MPlexHV    msPar;

  //============================================================================

  MkFinderFV() {}

  static constexpr int nnfv() { return NNFV; }
  static constexpr int nCands() { return ncands; }
  static constexpr int nSeeds() { return nseeds; }
  int operator()(int seed, int track) { return seed*ncands + track; }

  //----------------------------------------------------------------------------

  void InputTrack(const Track& track, int iseed, int offset, bool inputProp);
  void OutputTrack(std::vector<Track>& tracks, int itrack, int imp, bool outputProp) const;

  //----------------------------------------------------------------------------

  void SelectHitIndices(const LayerOfHits &layer_of_hits);

  void FindCandidates(const LayerOfHits &layer_of_hits, const FindingFoos &fnd_foos);

  void UpdateWithLastHit(const LayerOfHits &layer_of_hits, const FindingFoos &fnd_foos);

  void SelectBestCandidates(const LayerOfHits &layer_of_hits);
  int BestCandidate(int offset) const;

  void CopyOutParErr(std::vector<CombCandidate>& seed_cand_vec,
                     int N_proc, bool outputProp) const;

  //----------------------------------------------------------------------------

private:

  int XHitMax() const {
    int m = 0;
    for (int i = 0; i < NNFV; ++i) m = std::max(m, XHitSize[i]);
    return m;
  }

  int XHitMax(int offset) const {
    int m = 0;
    int j = offset*ncands;
    for (int i = 0; i < ncands; ++i, ++j) m = std::max(m, XHitSize[j]);
    return m;
  }

  void copy_in(const Track& trk, const int mslot, const int tslot)
  {
    Err[tslot].CopyIn(mslot, trk.errors().Array());
    Par[tslot].CopyIn(mslot, trk.parameters().Array());

    Chg  (mslot, 0, 0) = trk.charge();
    Chi2 (mslot, 0, 0) = trk.chi2();
    Label(mslot, 0, 0) = trk.label();

    NHits     (mslot, 0, 0) = trk.nTotalHits();
    NFoundHits(mslot, 0, 0) = trk.nFoundHits();
    std::copy(trk.BeginHitsOnTrack(), trk.EndHitsOnTrack(), HoTArrs[mslot]); 
  }

  void copy_out(Track& trk, const int mslot, const int tslot) const
  {
    Err[tslot].CopyOut(mslot, trk.errors_nc().Array());
    Par[tslot].CopyOut(mslot, trk.parameters_nc().Array());

    trk.setCharge(Chg  (mslot, 0, 0));
    trk.setChi2  (Chi2 (mslot, 0, 0));
    trk.setLabel (Label(mslot, 0, 0));

    trk.setNTotalHits(NHits     (mslot, 0, 0));
    trk.setNFoundHits(NFoundHits(mslot, 0, 0));
    std::copy(HoTArrs[mslot], & HoTArrs[mslot][NHits(mslot, 0, 0)], trk.BeginHitsOnTrack_nc());
  }

  void add_hit(const int mslot, int index, int layer)
  {
    int &n_tot_hits = NHits(mslot, 0, 0);
    int &n_fnd_hits = NFoundHits(mslot, 0, 0);

    if (n_tot_hits < Config::nMaxTrkHits)
    {
      HoTArrs[mslot][n_tot_hits++] = { index, layer };
      if (index >= 0) { ++n_fnd_hits; }
    }
    else
    {
      // printf("WARNING MkFinder::add_hit hit-on-track limit reached for label=%d\n", label_);

      const int LH = Config::nMaxTrkHits - 1;

      if (index >= 0)
      {
        if (HoTArrs[mslot][LH].index < 0)
          ++n_fnd_hits;
        HoTArrs[mslot][LH] = { index, layer };
      }
      else if (index == -2)
      {
        if (HoTArrs[mslot][LH].index >= 0)
          --n_fnd_hits;
        HoTArrs[mslot][LH] = { index, layer };
      }
    }
  }

  int num_invalid_hits(const int mslot) const
  {
    return NHits(mslot,0,0) - NFoundHits(mslot,0,0);
  }
};

using MkFinderFv = MkFinderFV<NN/8, 8>;
MkFinderFv* makeMkFinderFv();

#endif
