#ifndef CCLUSTERING_GMSTCLUSTERING_H
#define CCLUSTERING_GMSTCLUSTERING_H

#include <algorithm>
#include <memory>

#include <G4ExcitationHandler.hh>
#include <G4FermiPhaseDecay.hh>
#include <G4ReactionProductVector.hh>
#include <G4NucleiProperties.hh>

#include "MSTClustering.hh"

class CoordinateMSTClustering : public MSTClustering {
  public:
    CoordinateMSTClustering() = delete;
    CoordinateMSTClustering(bool consider_rep, bool extra_momentum, int stat_exen_type) : _consider_rep(consider_rep),
      _extra_momentum(extra_momentum), _stat_exen_type(stat_exen_type) {};
  private:

    static constexpr double nucleonAverMass = 0.93891875434*CLHEP::GeV;
    bool _consider_rep;
    bool _extra_momentum;
    uint32_t _stat_exen_type;

    std::vector<Edge> get_edges(const cola::EventData&) final;
    std::unique_ptr<cola::EventData> get_clusters(std::unique_ptr<cola::EventData>&&) final;
    cola::EventParticles _process_side(const cola::EventData&, cola::ParticleClass);

    G4FermiPhaseDecay phaseSpaceDecay;
};

#endif //CCLUSTERING_GMSTCLUSTERING_H