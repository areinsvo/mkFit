#ifndef _track_
#define _track_

#include "Hit.h"
#include "Matrix.h"
#include <vector>

struct TrackState
{
public:
  TrackState() : valid(true) {}
  SVector6 parameters;
  SMatrixSym66 errors;
  int charge;
  bool valid;
};

class Track
{
public:
  Track() {}

  Track(TrackState state, HitVec hits, float chi2) : state_(state), hits_(hits), chi2_(chi2) {}
  Track(int charge, SVector3 position, SVector3 momentum, SMatrixSym66 errors, HitVec hits, float chi2) {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = SVector6(position.At(0),position.At(1),position.At(2),momentum.At(0),momentum.At(1),momentum.At(2));
    state_.valid = true;
    hits_=hits;
    chi2_=chi2;
  }
  Track(int charge, SVector3 position, SVector3 momentum, SMatrixSym66 errors, HitVec hits, float chi2, HitVec initHits) {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = SVector6(position.At(0),position.At(1),position.At(2),momentum.At(0),momentum.At(1),momentum.At(2));
    state_.valid = true;
    hits_=hits;
    initHits_=initHits;
    chi2_=chi2;
  }
  Track(int charge, SVector6& parameters, SMatrixSym66& errors,HitVec hits, float chi2) {
    state_.charge=charge;
    state_.errors=errors;
    state_.parameters = parameters;
    state_.valid = true;
    hits_=hits;
    chi2_=chi2;
  }

  ~Track(){}

  int           charge() const {return state_.charge;}
  SVector3      position() const {return SVector3(state_.parameters[0],state_.parameters[1],state_.parameters[2]);}
  SVector3      momentum() const {return SVector3(state_.parameters[3],state_.parameters[4],state_.parameters[5]);}
  const SVector6&     parameters() const {return state_.parameters;}
  const SMatrixSym66& errors() const {return state_.errors;}
  const TrackState&   state() const {return state_;}
  float         chi2() const {return chi2_;}

  const HitVec& hitsVector() const {return hits_;}
  const HitVec& initHitsVector() const {return initHits_;}
  void addHit(Hit hit,float chi2) {hits_.push_back(hit);chi2_+=chi2;}
  void addInitHit(Hit hit);
  void resetHits() {hits_.clear();}
  unsigned int nHits() const {return hits_.size();}
  unsigned int SimTrackID() const;
  std::vector<unsigned int> SimTrackIDs() const;

  Track clone() const {return Track(state_,hits_,chi2_);}

private:
  TrackState state_;
  HitVec hits_;
  HitVec initHits_;
  float chi2_;

};

typedef std::vector<Track> TrackVec;
unsigned int getPhiPartition(float phi);
#endif
