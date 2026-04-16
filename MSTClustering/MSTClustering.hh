//
// Created by _amp_ on 10/21/24.
//

#ifndef CCLUSTERING_MSTCLUSTERING_H
#define CCLUSTERING_MSTCLUSTERING_H

#include <optional>
#include <memory>
#include <algorithm>
#include <vector>
#include <map>

#include "COLA.hh"

class MSTClustering : public cola::VConverter {
public:
    MSTClustering() = default;
    ~MSTClustering() override = default;
    MSTClustering& operator=(const MSTClustering&) = delete;
    MSTClustering& operator=(MSTClustering&&) = delete;
    MSTClustering(const MSTClustering&) = delete;
    MSTClustering(MSTClustering&&) = delete;

    std::unique_ptr<cola::EventData> operator()(std::unique_ptr<cola::EventData>&& data) final;

protected:
    using nPair = std::pair<cola::Particle*, cola::Particle*>;        // vertices pair
    struct Edge {                               // edge
        nPair vert;
        double size;
        cola::ParticleClass pClass;
        Edge(nPair vert_, double size_, cola::ParticleClass pClass_) : vert(vert_), size(size_), pClass(pClass_) {};
    };

    // a single Node with children.
    struct Node {
        explicit Node(cola::Particle* vertex) : height(0.), vertices(1, vertex) {}
        Node() : Node(0) {}
        Node(Node* first, Node* second, double height_) : height(height_),
        vertices(first->vertices), children(std::make_pair(first, second)) {
            vertices.insert(vertices.end(), second->vertices.begin(), second->vertices.end()); // append second vector
        }

        double height;
        std::vector<cola::Particle*> vertices;
        std::optional<std::pair<Node*, Node*>> children;
    };

    Node* rootA;
    Node* rootB;

    // it is reasonable to use the iterators for quick access to spectators even at this abstract level, so this class has it
    cola::EventParticles::iterator spectIterA;
    cola::EventParticles::iterator spectIterB;
    cola::EventParticles::iterator endIter;

private:

    void construct_trees(std::vector<Edge>&& edgeData, std::vector<std::unique_ptr<Node>>& nodes);

    // get full graph data from EventData
    virtual std::vector<Edge> get_edges(const cola::EventData&) = 0;
    // construct clusters
    virtual std::unique_ptr<cola::EventData> get_clusters(std::unique_ptr<cola::EventData>&&) = 0;
};



#endif // CCLUSTERING_MSTCLUSTERING_H
