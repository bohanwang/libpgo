/*
author: Bohan Wang
copyright to USC
*/

#include "triangleMeshSelfContactDetection.h"

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>

#include <vector>
#include <cmath>
#include <cfloat>

using namespace pgo;
using namespace pgo::Contact;
using namespace pgo::Mesh;

namespace ES = pgo::EigenSupport;

TriangleMeshSelfContactDetection::TriangleMeshSelfContactDetection(TriMeshRef mesh)
    : triangleMeshRef(mesh), triMeshNeighbor(triangleMeshRef) {
    bvTree.buildByInertiaPartition(triangleMeshRef);

    lastAABBs.resize(bvTree.getNumNodes());
    curAABBs.resize(bvTree.getNumNodes());
    for (int i = 0; i < bvTree.getNumNodes(); i++) {
        lastAABBs[i] = bvTree.getNodeBoundingBox(i);
    }

    for (int i = 0; i < triangleMeshRef.numTriangles(); i++) {
        triangles.push_back(triangleMeshRef.tri(i)[0]);
        triangles.push_back(triangleMeshRef.tri(i)[1]);
        triangles.push_back(triangleMeshRef.tri(i)[2]);
    }

    lastTriangleAABBs = triangleMeshRef.getTriangleBoundingBoxes();
    curTriangleAABBs  = lastTriangleAABBs;

    travseralFrontier0.reserve(100000);
    travseralFrontier1.reserve(100000);
}

void TriangleMeshSelfContactDetection::execute(const double* positions0, const double* positions1) {
    execute(positions0, positions1, 0.0);
}

void TriangleMeshSelfContactDetection::execute(const double* positions0, double activationDistance) {
    execute(positions0, nullptr, activationDistance);
}

void TriangleMeshSelfContactDetection::execute(const double* positions0, const double* positions1,
                                               double activationDistance) {
    TriMeshRef meshRef0(triangleMeshRef.numVertices(), positions0, triangleMeshRef.numTriangles(), triangles.data());
    bvTree.updateBoundingVolumes(meshRef0);

    for (int i = 0; i < bvTree.getNumNodes(); i++) {
        lastAABBs[i] = bvTree.getNodeBoundingBox(i);
    }
    lastTriangleAABBs = meshRef0.getTriangleBoundingBoxes();

    int isCCD = 0;
    if (positions1) {
        TriMeshRef meshRef1(triangleMeshRef.numVertices(), positions1, triangleMeshRef.numTriangles(),
                            triangles.data());
        bvTree.updateBoundingVolumes(meshRef1);

        for (int i = 0; i < bvTree.getNumNodes(); i++) {
            curAABBs[i] = bvTree.getNodeBoundingBox(i);
        }
        curTriangleAABBs = meshRef1.getTriangleBoundingBoxes();
        isCCD = 1;
    } else {
        curTriangleAABBs = lastTriangleAABBs;
    }

    travseralFrontier0.clear();
    travseralFrontier1.clear();

    lastFrontier    = &travseralFrontier0;
    currentFrontier = &travseralFrontier1;
    lastFrontier->emplace_back(bvTree.getRootNodeID(), bvTree.getRootNodeID());

    for (auto it = trianglePairBufferTLS.begin(); it != trianglePairBufferTLS.end(); ++it)
        it->clear();

    int inc = 0;
    while (lastFrontier->size()) {
        for (auto it = bvttNodeBufferTLS.begin(); it != bvttNodeBufferTLS.end(); ++it)
            it->clear();

        if (lastFrontier->size() < parallel_threshold) {
            for (const auto& node : *lastFrontier) {
                addressSingleNode(node, isCCD, activationDistance);
            }  // end for each node in last frontier
        }  // end single thread
        else {
            size_t count = (lastFrontier->size() + parallel_threshold - 1) / parallel_threshold;
            tbb::parallel_for((size_t)0, count, [&](size_t ci) {
                size_t start = ci * parallel_threshold;
                size_t end   = std::min(start + parallel_threshold, lastFrontier->size());
                for (size_t ni = start; ni < end; ni++)
                    addressSingleNode(lastFrontier->at(ni), isCCD, activationDistance);
            });
        }

        currentFrontier->clear();
        for (auto it = bvttNodeBufferTLS.begin(); it != bvttNodeBufferTLS.end(); ++it) {
            currentFrontier->insert(currentFrontier->end(), it->begin(), it->end());
        }

        std::swap(lastFrontier, currentFrontier);
        // std::cout << inc << ':' << lastFrontier->size() << std::endl;
        inc++;
        // std::cin.get();
    }

    candidateTrianglePairs.clear();
    for (auto it = trianglePairBufferTLS.begin(); it != trianglePairBufferTLS.end(); ++it) {
        candidateTrianglePairs.insert(candidateTrianglePairs.end(), it->begin(), it->end());
    }

    // address triangle pair ccd
    tbb::parallel_for((size_t)0, candidateTrianglePairs.size(),
                      [&](size_t ti) {
                          if (candidateTrianglePairs[ti].first > candidateTrianglePairs[ti].second)
                              std::swap(candidateTrianglePairs[ti].first, candidateTrianglePairs[ti].second);
                      },
                      tbb::static_partitioner());

    auto cmpFunc = [](const std::pair<int, int>& pa, const std::pair<int, int>& pb) -> bool {
        return std::less<std::pair<int, int>>()(pa, pb);
    };

    auto equalFunc = [](const std::pair<int, int>& pa, const std::pair<int, int>& pb) -> bool {
        return std::equal_to<std::pair<int, int>>()(pa, pb);
    };

    if (candidateTrianglePairs.size() > 1e5)
        tbb::parallel_sort(candidateTrianglePairs.begin(), candidateTrianglePairs.end(), cmpFunc);
    else
        std::sort(candidateTrianglePairs.begin(), candidateTrianglePairs.end(), cmpFunc);

    auto cdTriIt = std::unique(candidateTrianglePairs.begin(), candidateTrianglePairs.end(), equalFunc);
    candidateTrianglePairs.erase(cdTriIt, candidateTrianglePairs.end());

    // std::cout << candidateTrianglePairs.size() << std::endl;
}

void TriangleMeshSelfContactDetection::addressSingleNode(const BVTTNode& node, int isCCD, double activationDistance) {
    BVTTNodeBuffer& nodeBuffer = bvttNodeBufferTLS.local();
    // if they are the same node
    if (node.nodeA == node.nodeB) {
        // if they are not leaves
        if (bvTree.getNumElements(node.nodeA) == 0) {
            expandBothNodes(node.nodeA, node.nodeB, isCCD, activationDistance, nodeBuffer);
        }
        // else if they are leaves
        else {
            assignTriangles(node.nodeA, node.nodeB, isCCD, activationDistance, trianglePairBufferTLS.local());
        }
    }
    // else if they are different nodes, we have three different case
    else {
        // if they are not leaves
        if (bvTree.getNumElements(node.nodeA) == 0 && bvTree.getNumElements(node.nodeB) == 0) {
            expandBothNodes(node.nodeA, node.nodeB, isCCD, activationDistance, nodeBuffer);
        }
        // else if A is not leaf but B is leaf
        else if (bvTree.getNumElements(node.nodeA) == 0 && bvTree.getNumElements(node.nodeB) != 0) {
            expandFirstNodes(node.nodeA, node.nodeB, isCCD, activationDistance, nodeBuffer);
        }
        // else if B is not leaf but A is leaf
        else if (bvTree.getNumElements(node.nodeA) != 0 && bvTree.getNumElements(node.nodeB) == 0) {
            expandFirstNodes(node.nodeB, node.nodeA, isCCD, activationDistance, nodeBuffer);
        }
        // else both are leaf node
        else {
            assignTriangles(node.nodeA, node.nodeB, isCCD, activationDistance, trianglePairBufferTLS.local());
        }
    }  // end else nodeA==nodeB
}

void TriangleMeshSelfContactDetection::assignTriangles(int nodeA, int nodeB, int isCCD, double activationDistance,
                                                       TrianglePairBuffer& buf) {
    for (int i = 0; i < bvTree.getNumElements(nodeA); i++) {
        int triIdA = bvTree.getElementID(nodeA, i);
        // Vec3i neighborA = triMeshNeighbor.getTriangleNeighbors(triIdA);

        for (int j = (nodeA == nodeB ? i : 0); j < bvTree.getNumElements(nodeB); j++) {
            int triIdB = bvTree.getElementID(nodeB, j);
            // if B is not a neighbor of A
            if (triIdA != triIdB &&            // not same triangle
                !shareVertex(triIdA, triIdB) &&
                trianglesWithinQueryBand(triIdA, triIdB, isCCD, activationDistance))  // not vertex/edge neighbor
            {
                buf.emplace_back(triIdA, triIdB);
            }
        }
    }
}

void TriangleMeshSelfContactDetection::expandBothNodes(int nodeA, int nodeB, int isCCD, double activationDistance,
                                                       BVTTNodeBuffer& buf) {
    // for each pair of children nodes
    for (int i = 0; i < bvTree.getNumChildNodes(nodeA); i++) {
        int childA = bvTree.getChildNodeID(nodeA, i);

        for (int j = (nodeA == nodeB ? i : 0); j < bvTree.getNumChildNodes(nodeB); j++) {
            int childB = bvTree.getChildNodeID(nodeB, j);

            // if they are the same node, they are considered intersected
            if (childA == childB)
                buf.emplace_back(childA, childB);
            else if (nodesWithinQueryBand(childA, childB, isCCD, activationDistance)) {
                buf.emplace_back(childA, childB);
            }
        }
    }
}

void TriangleMeshSelfContactDetection::expandFirstNodes(int nodeA, int nodeB, int isCCD, double activationDistance,
                                                        BVTTNodeBuffer& buf) {
    // for each pair of children nodes
    for (int i = 0; i < bvTree.getNumChildNodes(nodeA); i++) {
        int childA = bvTree.getChildNodeID(nodeA, i);

        // if they are the same node, they are considered intersected
        if (childA == nodeB)
            buf.emplace_back(childA, nodeB);
        else if (nodesWithinQueryBand(childA, nodeB, isCCD, activationDistance)) {
            buf.emplace_back(childA, nodeB);
        }
    }
}

bool TriangleMeshSelfContactDetection::nodesWithinQueryBand(int nodeA, int nodeB, int isCCD,
                                                            double activationDistance) const {
    const double band = std::max(activationDistance, 0.0);

    if (band <= 0.0) {
        if (isCCD) {
            return CCDIntersect(lastAABBs[nodeA], curAABBs[nodeA], lastAABBs[nodeB], curAABBs[nodeB]);
        }

        return DCDIntersect(lastAABBs[nodeA], lastAABBs[nodeB]);
    }

    if (isCCD) {
        return AABBWithinBand(SweptAABB(lastAABBs[nodeA], curAABBs[nodeA]),
                              SweptAABB(lastAABBs[nodeB], curAABBs[nodeB]), band);
    }

    return AABBWithinBand(lastAABBs[nodeA], lastAABBs[nodeB], band);
}

bool TriangleMeshSelfContactDetection::trianglesWithinQueryBand(int triA, int triB, int isCCD,
                                                                double activationDistance) const {
    const double band = std::max(activationDistance, 0.0);

    if (band <= 0.0) {
        if (isCCD) {
            return CCDIntersect(lastTriangleAABBs[triA], curTriangleAABBs[triA], lastTriangleAABBs[triB],
                                curTriangleAABBs[triB]);
        }

        return DCDIntersect(lastTriangleAABBs[triA], lastTriangleAABBs[triB]);
    }

    if (isCCD) {
        return AABBWithinBand(SweptAABB(lastTriangleAABBs[triA], curTriangleAABBs[triA]),
                              SweptAABB(lastTriangleAABBs[triB], curTriangleAABBs[triB]), band);
    }

    return AABBWithinBand(lastTriangleAABBs[triA], lastTriangleAABBs[triB], band);
}

bool TriangleMeshSelfContactDetection::shareVertex(int triA, int triB) {
    Vec3i triAVtxId = triangleMeshRef.tri(triA);
    Vec3i triBVtxId = triangleMeshRef.tri(triB);

    return triAVtxId[0] == triBVtxId[0] || triAVtxId[0] == triBVtxId[1] || triAVtxId[0] == triBVtxId[2] ||
           triAVtxId[1] == triBVtxId[0] || triAVtxId[1] == triBVtxId[1] || triAVtxId[1] == triBVtxId[2] ||
           triAVtxId[2] == triBVtxId[0] || triAVtxId[2] == triBVtxId[1] || triAVtxId[2] == triBVtxId[2];
}

bool TriangleMeshSelfContactDetection::CCDIntersect(const BoundingBox& AABB1Start, const BoundingBox& AABB1End,
                                                    const BoundingBox& AABB2Start, const BoundingBox& AABB2End) {
    double collLo[3], collHi[3];

    for (int dim = 0; dim < 3; dim++) {
        double LO1S = AABB1Start.bmin()[dim];
        double HI1S = AABB1Start.bmax()[dim];
        double LO1E = AABB1End.bmin()[dim];
        double HI1E = AABB1End.bmax()[dim];

        double LO2S = AABB2Start.bmin()[dim];
        double HI2S = AABB2Start.bmax()[dim];
        ;
        double LO2E = AABB2End.bmin()[dim];
        double HI2E = AABB2End.bmax()[dim];
        ;

        double LOS = 0.0, HIS = 0.0, LOE = 0.0, HIE = 0.0;

        // check if collision at t=0
        // [LO1S, HI1S] , [LO2S, HI2S]
        int collisionAtZero = 0;
        if (LO2S > HI1S) {
            LOS = HI1S;
            LOE = HI1E;
            HIS = LO2S;
            HIE = LO2E;
        } else if (LO1S > HI2S) {
            LOS = HI2S;
            LOE = HI2E;
            HIS = LO1S;
            HIE = LO1E;
        } else {
            // collision at t = 0
            collLo[dim]     = 0.0;
            collisionAtZero = 1;
        }

        if (!collisionAtZero) {
            // determine intersection point of LOS, LOE to HIS, HIE
            // we know LOS < HIS
            if (LOE < HIE)
                collLo[dim] = DBL_MAX;
            // compute intersection point
            // (1-t) * LOS + t * LOE = (1-t) * HIS + t * HIE
            collLo[dim] = (HIS - LOS) / (LOE - LOS + HIS - HIE);
        }

        // check if collision at t=1
        // [LO1E, HI1E] , [LO2E, HI2E]
        int collisionAtOne = 0;
        if (LO2E > HI1E) {
            LOS = HI1S;
            LOE = HI1E;
            HIS = LO2S;
            HIE = LO2E;
        } else if (LO1E > HI2E) {
            LOS = HI2S;
            LOE = HI2E;
            HIS = LO1S;
            HIE = LO1E;
        } else {
            // collision at t = 1
            collHi[dim]    = 1.0;
            collisionAtOne = 1;
        }

        if (!collisionAtOne) {
            // determine intersection point of LOS, LOE to HIS, HIE
            // we know LOS < HIS
            if (LOS < HIS)
                collHi[dim] = -DBL_MAX;
            // compute intersection point
            // (1-t) * LOS + t * LOE = (1-t) * HIS + t * HIE
            collHi[dim] = (HIS - LOS) / (LOE - LOS + HIS - HIE);
        }
    }

    double lo = std::max(std::max(collLo[0], collLo[1]), collLo[2]);
    double hi = std::min(std::min(collHi[0], collHi[1]), collHi[2]);

    return (lo <= hi);
}

bool TriangleMeshSelfContactDetection::DCDIntersect(const BoundingBox& AABB1, const BoundingBox& AABB2) {
    return AABB1.intersect(AABB2);
}

bool TriangleMeshSelfContactDetection::AABBWithinBand(const BoundingBox& AABB1, const BoundingBox& AABB2, double band) {
    const double clampedBand = std::max(band, 0.0);

    for (int dim = 0; dim < 3; ++dim) {
        double separation = 0.0;
        if (AABB1.bmax()[dim] < AABB2.bmin()[dim]) {
            separation = AABB2.bmin()[dim] - AABB1.bmax()[dim];
        } else if (AABB2.bmax()[dim] < AABB1.bmin()[dim]) {
            separation = AABB1.bmin()[dim] - AABB2.bmax()[dim];
        }

        if (separation > clampedBand) {
            return false;
        }
    }

    return true;
}

BoundingBox TriangleMeshSelfContactDetection::SweptAABB(const BoundingBox& start, const BoundingBox& end) {
    BoundingBox swept = start;
    swept.expand(end);
    return swept;
}
