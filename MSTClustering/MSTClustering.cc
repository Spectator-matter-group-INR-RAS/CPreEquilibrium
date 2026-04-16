//
// Created by _amp_ on 10/21/24.
//

#include "MSTClustering.hh"
#include "G4NucleiProperties.hh"
#include "Repulsion.hh"
#include "G4ExcitationHandler.hh"
#include "ExcitationEnergy.hh"
#include <limits>
#include <stack>

void MSTClustering::construct_trees(std::vector<Edge> &&edgeData, std::vector<std::unique_ptr<Node>>& nodes)
{
    std::unordered_map<cola::Particle*, Node*> treeA;
    std::unordered_map<cola::Particle*, Node*> treeB;

    // initialize trees
    for (auto iter = spectIterA; iter < spectIterB; ++iter) {
        nodes.emplace_back(std::make_unique<Node>(&(*iter)));
        treeA.emplace(&(*iter), nodes.back().get());
    }
    for (auto iter = spectIterB; iter < endIter; ++iter) {
        nodes.emplace_back(std::make_unique<Node>(&(*iter)));
        treeB.emplace(&(*iter), nodes.back().get());
    }
    
    // sort the edges for hierarchical tree
    std::sort(edgeData.begin(), edgeData.end(), [](const Edge& l, const Edge& r) { return l.size < r.size; });
    // merge nodes into complete trees using edgeData
    for (const auto& it : edgeData)
    {
        auto v1 = it.vert.first;
        auto v2 = it.vert.second;
        switch (it.pClass)
        {
        case cola::ParticleClass::spectatorA:
            if (treeA[v1] != treeA[v2]) {
                nodes.emplace_back(std::make_unique<Node>(treeA[v1], treeA[v2], it.size));
                for (const auto vertex : nodes.back()->vertices) {
                    treeA[vertex] = nodes.back().get();
                }
            }
            break;
        case cola::ParticleClass::spectatorB:
            if (treeB[v1] != treeB[v2]) {
                nodes.emplace_back(std::make_unique<Node>(treeB[v1], treeB[v2], it.size));
                for (const auto vertex : nodes.back()->vertices) {
                    treeB[vertex] = nodes.back().get();
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
    std::sort(data->particles.begin(), data->particles.end(), [](cola::Particle l, cola::Particle r) {return l.pClass < r.pClass;});
    spectIterA = std::find_if(data->particles.begin(), data->particles.end(), [](cola::Particle p) {return p.pClass == cola::ParticleClass::spectatorA;});
    spectIterB = std::find_if(data->particles.begin(), data->particles.end(), [](cola::Particle p) {return p.pClass == cola::ParticleClass::spectatorB;});
    endIter = data->particles.end();
    // construct trees
    std::vector<std::unique_ptr<Node>> nodeHolder;
    construct_trees(get_edges(*data), nodeHolder);
    // divide trees and process resulting pre-fragments
    return get_clusters(std::move(data));
}
