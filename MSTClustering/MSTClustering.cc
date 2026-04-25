//
// Created by _amp_ on 10/21/24.
//

#include <G4NucleiProperties.hh>
#include <G4ExcitationHandler.hh>
#include <limits>
#include <stack>

#include "Repulsion.hh"
#include "ExcitationEnergy.hh"

#include "MSTClustering.hh"

void MSTClustering::construct_trees(std::vector<Edge>&& edgeData, std::vector<Node>& nodes)
{
    std::unordered_map<cola::Particle*, Node*> treeA;
    std::unordered_map<cola::Particle*, Node*> treeB;
    treeA.reserve(std::distance(spectIterA, spectIterB));
    treeB.reserve(std::distance(spectIterB, endIter));

    // no resizes are allowed (otherwise segfault)
    nodes.reserve(2 * (treeA.bucket_count() + treeB.bucket_count()));

    // initialize trees
    for (auto iter = spectIterA; iter < spectIterB; ++iter) {
        nodes.emplace_back(*iter);
        treeA.emplace(&(*iter), &nodes.back());
    }
    for (auto iter = spectIterB; iter < endIter; ++iter) {
        nodes.emplace_back(*iter);
        treeB.emplace(&(*iter), &nodes.back());
    }

    // sort the edges for hierarchical tree
    std::sort(edgeData.begin(), edgeData.end(), [](const Edge& l, const Edge& r) { return l.size < r.size; });

    // merge nodes into complete trees using edgeData
    for (const auto& edge : edgeData)
    {
        auto v1 = edge.vert.first;
        auto v2 = edge.vert.second;
        switch (edge.pClass)
        {
        case cola::ParticleClass::spectatorA:
            if (treeA[v1] != treeA[v2]) {
                nodes.emplace_back(treeA[v1], treeA[v2], edge.size);
                for (const auto vertex : nodes.back().vertices) {
                    treeA[vertex] = &nodes.back();
                }
            }
            break;
        case cola::ParticleClass::spectatorB:
            if (treeB[v1] != treeB[v2]) {
                nodes.emplace_back(treeB[v1], treeB[v2], edge.size);
                for (const auto vertex : nodes.back().vertices) {
                    treeB[vertex] = &nodes.back();
                }
            }
            break;
        default:
            throw(std::logic_error("An edge between non-Spectators have been formed!"));
        }
    }

    // set root nodes

    rootA = !treeA.empty() ? treeA.begin()->second : nullptr;
    rootB = !treeB.empty() ? treeB.begin()->second : nullptr;
}

std::unique_ptr<cola::EventData> MSTClustering::operator()(std::unique_ptr<cola::EventData>&& data) {
    // sort particles by pClass (spectators are last)
    std::sort(data->particles.begin(), data->particles.end(), [](const cola::Particle& l, const cola::Particle& r) { return l.pClass < r.pClass; });
    spectIterA = std::find_if(data->particles.begin(), data->particles.end(), [](cola::Particle p) {return p.pClass == cola::ParticleClass::spectatorA;});
    spectIterB = std::find_if(data->particles.begin(), data->particles.end(), [](cola::Particle p) {return p.pClass == cola::ParticleClass::spectatorB;});
    endIter = data->particles.end();
    // construct trees
    std::vector<Node> nodes; // stores all nodes in trees (improves cpu cache hits)
    construct_trees(get_edges(*data), nodes);
    // divide trees and process resulting pre-fragments
    return get_clusters(std::move(data));
}
