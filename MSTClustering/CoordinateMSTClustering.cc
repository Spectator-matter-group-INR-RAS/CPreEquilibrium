#include <limits>
#include <queue>
#include <G4SystemOfUnits.hh>

#include "Repulsion.hh"
#include "ExcitationEnergy.hh"

#include "CoordinateMSTClustering.hh"

// constants for implementation
constexpr double eps0 = 2.17 * MeV;
constexpr double alphaPow = -1.02;
constexpr double d0 = 2.7;
constexpr double aColRel = 5.315;
constexpr double aSurfRel = 0.017;
constexpr double aVol = 0.054;

constexpr double SpecAa = 0;
constexpr double SpecAb = 0;

constexpr double a_opt = 2.243;
constexpr double b_opt = 3.183 * MeV;
constexpr double c_opt = 0.99;
constexpr double d_opt = 0.29041;

std::unique_ptr<cola::EventData> CoordinateMSTClustering::get_clusters(std::unique_ptr<cola::EventData> &&data)
{
    cola::EventParticles clustersA;
    cola::EventParticles clustersB;

    // get clusters
    if (rootA != nullptr)
    {
        clustersA = _process_side(*data, cola::ParticleClass::spectatorA);
    }
    if (rootB != nullptr)
    {
        clustersB = _process_side(*data, cola::ParticleClass::spectatorB);
    }

    // erase spectator nucleons
    data->particles.erase(spectIterA != endIter ? spectIterA : spectIterB, endIter);

    // append clusters
    if (rootA != nullptr)
    {
        data->particles.insert(data->particles.end(), std::make_move_iterator(clustersA.begin()), std::make_move_iterator(clustersA.end()));
    }
    if (rootB != nullptr)
    {
        data->particles.insert(data->particles.end(), std::make_move_iterator(clustersB.begin()), std::make_move_iterator(clustersB.end()));
    }

    return std::move(data);
}

std::vector<MSTClustering::Edge> CoordinateMSTClustering::get_edges(const cola::EventData&)
// Notice that we don't use EventData in this implementation since we have the needed iterators. The possibility is still there though
{
    std::vector<Edge> edges;
    edges.reserve(std::pow(std::distance(spectIterA, spectIterB), 2) + std::pow(std::distance(spectIterB, endIter), 2));

    // particle vector is sorted, process spectatorA nucleons (check for no spectatorA nucleons)
    if (spectIterA != endIter)
    {
        for (auto iter = spectIterA; iter != spectIterB; ++iter)
        {
            for (auto jter = iter + 1; jter != spectIterB; ++jter)
            {
                auto delta = iter->position - jter->position;
                edges.emplace_back(std::make_pair(&(*iter), &(*jter)), std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z), cola::ParticleClass::spectatorA);
            }
        }
    }
    // repeat for spectatorB nucleons
    for (auto iter = spectIterB; iter != endIter; ++iter)
    {
        for (auto jter = iter + 1; jter != endIter; ++jter)
        {
            auto delta = iter->position - jter->position;
            edges.emplace_back(std::make_pair(&(*iter), &(*jter)), std::sqrt(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z), cola::ParticleClass::spectatorB);
        }
    }
    return edges;
}

double get_cd(double Ex, uint32_t A)
{
    const auto Afloat = static_cast<double>(A);
    if (Ex / Afloat < eps0)
    {
        return d0;
    }
    double ex = Ex / Afloat;
    double dep = std::exp(-std::pow((ex / b_opt), a_opt)) * c_opt + d_opt;
    return d0 * std::pow(dep, 1. / 3.);
}

cola::LorentzVector ToColaLorentzVector(const G4LorentzVector& lv)
{
    return {
        lv.e(),
        lv.x(),
        lv.y(),
        lv.z(),
    };
}

cola::EventParticles CoordinateMSTClustering::_process_side(const cola::EventData& data, cola::ParticleClass side)
{
    cola::EventParticles clusters;

    auto& root = side == cola::ParticleClass::spectatorA ? rootA : rootB;
    auto bIter = side == cola::ParticleClass::spectatorA ? spectIterA : spectIterB;
    auto eIter = side == cola::ParticleClass::spectatorA ? spectIterB : endIter;

    uint32_t sourceA = cola::pdgToAZ(side == cola::ParticleClass::spectatorA ? data.iniState.pdgCodeA : data.iniState.pdgCodeB).first;
    // boost to rest frame for each set of spectators
    cola::LorentzVector pNucleus = {0.0, 0.0, 0.0, 0.0};
    for (auto particleIt = bIter; particleIt != eIter; ++particleIt)
    {
        pNucleus += particleIt->momentum;
    }
    for (auto particleIt = bIter; particleIt != eIter; ++particleIt)
    {
        particleIt->momentum.boost(-pNucleus);
    }
    const auto count = std::distance(bIter, eIter);

    // get excitation energy
    double exEn = ExcitationEnergy(_stat_exen_type, sourceA).GetEnergy(count);

    // at this point the construct_tree() method has already built up MST trees for both sides

    double cd = get_cd(exEn, count);
    auto unprocessed = std::queue<Node*>();
    unprocessed.push(root);

    while (!unprocessed.empty()) {
        auto* topView = unprocessed.front();
        if (topView->height <= cd)
        {
            cola::AZ clusterAZ = {0, 0};
            cola::LorentzVector position{0,0,0,0}, momentum{0,0,0,0};

            for (const auto* nucleon : topView->vertices)
            {
                cola::AZ componentAZ = nucleon->getAZ();
                clusterAZ.first += componentAZ.first;
                clusterAZ.second += componentAZ.second;
                position += nucleon->position;
                momentum += nucleon->momentum;
            }
            cola::Particle cluster;
            cluster.position = position / topView->vertices.size();
            cluster.momentum = momentum;
            cluster.pdgCode = cola::AZToPdg(clusterAZ);
            cluster.pClass = side;
            clusters.push_back(cluster);
        } else if (topView->children.has_value())
        {
            unprocessed.push(topView->children.value().first);
            unprocessed.push(topView->children.value().second);
        }
        unprocessed.pop();
    }

    // now that we have defined clusters, we need to calculate additional momentum

    if (_extra_momentum && clusters.size() > 1) {
        std::vector<double> mass;
        mass.reserve(clusters.size());
        double totalMass = .0;

        for (const auto& cluster: clusters)
        {
            mass.push_back(G4NucleiProperties::GetNuclearMass(static_cast<G4int>(cluster.getAZ().first), static_cast<G4int>(cluster.getAZ().second)) + exEn * cluster.getAZ().first / sourceA);
            totalMass += mass.back();
        }

        totalMass += 1e-5*MeV; // fix for segfault
        const auto extra_momentum = phaseSpaceDecay.CalculateDecay(G4LorentzVector(totalMass, {}), mass);

        for (size_t i = 0; i < clusters.size(); ++i)
        {
            clusters[i].momentum += ToColaLorentzVector(extra_momentum.at(i));
        }
    }

    if (_consider_rep)
    {
        RepulsionStage::CalculateRepulsion(std::move(clusters));
    }

    // make mass consistent with G4 tables and boost back from rest frame
    for (auto& cluster: clusters)
    {
        cluster.momentum.e = std::sqrt(
            std::pow(G4NucleiProperties::GetNuclearMass(static_cast<G4int>(cluster.getAZ().first), static_cast<G4int>(cluster.getAZ().second)), 2) +
            cluster.momentum.spatialPart().mag2());
        cluster.momentum.boost(pNucleus);
    }

    return clusters;
}

