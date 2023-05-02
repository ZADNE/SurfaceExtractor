/*!
 *  @author Burkalo Boris (xburka00)
 */

/*
* Code in this file is inspired by code in this repository:
* https://github.com/sandialabs/miniIsosurface/tree/master/flyingEdges
*
* Algorithm in this file is based on method from Schroeder et al. introduced
* in the following article:
* https://ieeexplore.ieee.org/document/7348069
*/

#include <SurfaceExtractor/FlyingEdges.hpp>
#include <SurfaceExtractor/Volume.hpp>
#include <SurfaceExtractor/FlyingEdgesCases.hpp>

#include <iostream>
#include <cmath>
#include <limits>
#include "FlyingEdges.hpp"

namespace SurfaceExtractor
{

float FlyingEdges::m_isovalue = 31.f / 255.f;
float FlyingEdges::m_spacing = 1.f;

void FlyingEdges::extractSurface(const Volume& vol,
    std::vector<glm::vec3>& vertices, std::vector<glm::uvec3>& triangles, float isoValue)
{
    m_isovalue = isoValue;

    // Pass 1
    std::vector<uint8_t> edgeCases((vol.dims.x - 1) * vol.dims.y * vol.dims.z, 0);
    std::vector<Edge> edges(vol.dims.y*vol.dims.z);

    pass1(vol, edgeCases, edges);
    //

    // Pass 2
    std::vector<uint8_t> cubeCases((vol.dims.x - 1) * (vol.dims.y - 1) * (vol.dims.z - 1), 0);
    std::vector<unsigned> triangleCount((vol.dims.y - 1) * (vol.dims.z - 1), 0);

    pass2(vol, edgeCases, edges, cubeCases, triangleCount);
    //

    // Pass 3
    unsigned verticesTot = 0;
    unsigned trianglesTot = 0;

    pass3(vol, edges, triangleCount, trianglesTot, verticesTot);
    //

    // Pass 4
    triangles.resize(trianglesTot);
    vertices.resize(verticesTot);

    pass4(vol, edges, triangleCount, cubeCases, vertices, triangles);
    //

    std::cout << verticesTot << std::endl;
    std::cout << trianglesTot << std::endl;
}

void FlyingEdges::pass1(const Volume& volume, std::vector<uint8_t>& edgeCases, std::vector<Edge>& edges)
{
    // Build the edge case vector.
    for (size_t k = 0; k < volume.dims.z; k++)
    {
        for (size_t j = 0; j < volume.dims.y; j++)
        {
            for (size_t i = 0; i < volume.dims.x - 1; i++)
            {
                size_t curId = k * volume.dims.y * (volume.dims.x - 1) + j * (volume.dims.x - 1) + i;

                uint8_t cur = ((float)volume.at(i, j, k) / 255.f) < m_isovalue;
                uint8_t next = ((float)volume.at(i + 1, j, k) / 255.f) < m_isovalue;
                
                // Define the case of the edge:
                // 0 if both of the vertices on the edge are greater than ISOVALUE
                // 1 if the next is greater
                // 2 if current is greater
                // 3 if neither of the vertices is greater 
                edgeCases[curId] = (next << 1) | (cur);
            }
        }
    }

    // Find limits for each edge.
    // These two loops are separated because the limits loop
    // needs info about other edge cases.
    for (size_t k = 0; k < volume.dims.z; k++)
    {
        for (size_t j = 0; j < volume.dims.y; j++)
        {
            int xL = volume.dims.x;
            int xR = 0;

            for (size_t i = 0; i < volume.dims.x - 1; i++)
            {
                bool updateLimits = false;

                size_t curId = k * volume.dims.y * (volume.dims.x - 1) + j * (volume.dims.x - 1) + i;
                uint8_t curEdgeCase = edgeCases[curId];
                if (curEdgeCase == 1 || curEdgeCase == 2)
                {
                    updateLimits = true;
                }
                else 
                {
                    // Checking if either of the a_ijk is cut (if x, y or z of E is cut):
                    //   |  /
                    // y | / z
                    //   |/
                    //   +-------+------------>
                    //   |   x   |       E
                    if (j != volume.dims.y - 1)
                    {
                        size_t curIdY = k * volume.dims.y * (volume.dims.x - 1) + (j + 1) * (volume.dims.x - 1) + i;
                        uint8_t curEdgeCaseY = edgeCases[curIdY];
                        
                        // Check the cases of each of the points on x, y, z
                        if ((curEdgeCase + curEdgeCaseY) % 2 == 1)
                        {
                            updateLimits = true;
                        }
                    }
                    else if (k != volume.dims.z - 1)
                    {
                        size_t curIdZ = (k + 1) * volume.dims.y * (volume.dims.x - 1) + j * (volume.dims.x - 1) + i;
                        uint8_t curEdgeCaseZ = edgeCases[curIdZ];
                        
                        if ((curEdgeCase + curEdgeCaseZ) % 2 == 1)
                        {
                            updateLimits = true;
                        }
                    }
                }
                
                if (updateLimits)
                {
                    if (xL == volume.dims.x)
                        xL = i;
                    
                    xR = i + 1;
                }
            }

            edges[k * volume.dims.y + j].xL = xL;
            edges[k * volume.dims.y + j].xR = xR;
        }
    }
}

void FlyingEdges::pass2(const Volume& volume, std::vector<uint8_t>& edgeCases,
        std::vector<Edge>& edges, std::vector<uint8_t>& cubeCases, std::vector<unsigned>& triangleCount)
{
    bool found = false;
    for (size_t k = 0; k < volume.dims.z - 1; k++)
    {
        for (size_t j = 0; j < volume.dims.y - 1; j++)
        {
            // Current edge of the cube.
            Edge& currentEdge = edges.at(k * volume.dims.y + j);

            // Other three X edges that are on the cube.
            Edge& otherEdge1 = edges.at(k * volume.dims.y + (j + 1));
            Edge& otherEdge2 = edges.at((k + 1) * volume.dims.y + j);
            Edge& otherEdge3 = edges.at((k + 1) * volume.dims.y + (j + 1));

            int xL = std::min(currentEdge.xL, std::min(otherEdge1.xL, std::min(otherEdge2.xL, otherEdge3.xL)));
            int xR = std::max(currentEdge.xR, std::max(otherEdge1.xR, std::max(otherEdge2.xR, otherEdge3.xR)));

            if (xR < xL)
                continue;

            int yEnd = j == (volume.dims.y - 2);
            int zEnd = k == (volume.dims.z - 2);

            for (size_t i = xL; i < xR; i++)
            {
                int xEnd = i == (volume.dims.x - 2);

                // assign marching cube cases
                std::array<uint8_t, 4> cubeEdgeCases = 
                {
                    edgeCases.at(k * volume.dims.y * (volume.dims.x - 1) + j * (volume.dims.x - 1) + i),
                    edgeCases.at(k * volume.dims.y * (volume.dims.x - 1) + (j + 1) * (volume.dims.x - 1) + i),
                    edgeCases.at((k + 1) * volume.dims.y * (volume.dims.x - 1) + j * (volume.dims.x - 1) + i),
                    edgeCases.at((k + 1) * volume.dims.y * (volume.dims.x - 1) + (j + 1) * (volume.dims.x - 1) + i)
                };

                uint8_t cubeCase = getCubeCase(cubeEdgeCases);

                cubeCases[k * (volume.dims.y - 1) * (volume.dims.x - 1) + j * (volume.dims.x - 1) + i] = cubeCase;

                if (cubeCase == 0 || cubeCase == 255)
                    continue;
                
                int gridEdgeCondition = xEnd | (yEnd << 1) | (zEnd << 2);

                std::vector<EvalEdge> evalEdges = 
                {
                    EvalEdge(0, currentEdge.xIntersects,  0, 0,  false),
                    EvalEdge(0, currentEdge.yIntersects,  0, 3,  false),
                    EvalEdge(0, currentEdge.zIntersects,  0, 8,  false),
                    EvalEdge(0, currentEdge.yIntersects,  1, 1,  false),
                    EvalEdge(0, currentEdge.zIntersects,  1, 9,  false),
                    EvalEdge(0, otherEdge1.xIntersects, 2, 2,  false),
                    EvalEdge(0, otherEdge1.zIntersects, 2, 10, false),
                    EvalEdge(0, otherEdge2.xIntersects,  4, 4,  false),
                    EvalEdge(0, otherEdge2.yIntersects,  4, 7,  false),
                    EvalEdge(0, otherEdge1.zIntersects, 3, 11, false),
                    EvalEdge(0, otherEdge2.yIntersects,  5, 5,  false),
                    EvalEdge(0, otherEdge3.xIntersects,  6, 6,  false)
                };

                // Count number of triangles for this row of cubes based on the table value
                // corresponding to the cube case.
                triangleCount[k * (volume.dims.y - 1) + j] += config::triangleCnt[cubeCase];

                // Count the number of intersections in x, y and z (namely 0, 3 and 8).
                // Other are counted in only when the cube is on the edge of the voxel 
                // grid.
                // the cube is oriented in the isCut like:
                //          6
                //      +--------+
                //   10/|     11/|
                //    / | 2    / |
                //   +--------+  | 5
                //   | 7|  4  |1 |
                // 3 |  +-----|--+
                //   |8/      | /
                //   |/       |/ 9
                //   +--------+
                //       0
                // X edges = [0, 2,  4,  6]
                // Y edges = [1, 3,  5,  7]
                // Z edges = [8, 9, 10, 11]
                // CE = [0, 3, 8]
                // OE1 = [2, 10, 11]

                for (auto ee : evalEdges)
                {
                    if ((gridEdgeCondition & ee.condition) == ee.condition)
                        ee.extraCnt += config::cubeEdgeCuts[cubeCase][ee.id];
                }
            }
        }
    }
}

void FlyingEdges::pass3(const Volume& volume, std::vector<Edge>& edges,
    std::vector<unsigned>& triangleCount, unsigned& trianglesTot, unsigned& verticesTot)
{
    // Mark the indices of triangles and vertices for each edge.
    size_t trianglesTotal = 0;
    size_t verticesTotal = 0;
    for (size_t k = 0; k < volume.dims.z; k++)
    {
        for (size_t j = 0; j < volume.dims.y; j++)
        {
            if (k < volume.dims.z - 1 && j < volume.dims.y - 1)
            {
                size_t prevSum = trianglesTotal;
                unsigned edgeTriCount = triangleCount[k * (volume.dims.y - 1) + j];

                triangleCount[k * (volume.dims.y - 1) + j] = trianglesTotal;
                trianglesTotal += edgeTriCount;
            }

            size_t prevSum = edges[k * volume.dims.y + j].xIntersects;
            edges[k * volume.dims.y + j].xIntersects = verticesTotal;
            verticesTotal += prevSum;

            prevSum = edges[k * volume.dims.y + j].yIntersects;
            edges[k * volume.dims.y + j].yIntersects = verticesTotal;
            verticesTotal += prevSum;

            prevSum = edges[k * volume.dims.y + j].zIntersects;
            edges[k * volume.dims.y + j].zIntersects = verticesTotal;
            verticesTotal += prevSum;

        }
    }

    trianglesTot = trianglesTotal;
    verticesTot = verticesTotal;
}

void FlyingEdges::pass4(const Volume& volume, const std::vector<Edge>& edges, const std::vector<unsigned>& triangleCount,
    const std::vector<uint8_t>& cubeCases, std::vector<glm::vec3>& vertices, std::vector<glm::uvec3>& triangles)
{
    // Build the resulting indices and vertices.
    for (size_t k = 0; k < volume.dims.z - 1; k++)
    {
        for (size_t j = 0; j < volume.dims.y - 1; j++)
        {
            const Edge& currentEdge = edges.at(k * volume.dims.y + j);
            const Edge& otherEdge1  = edges.at(k * volume.dims.y + (j + 1));
            const Edge& otherEdge2  = edges.at((k + 1) * volume.dims.y + j);
            const Edge& otherEdge3  = edges.at((k + 1) * volume.dims.y + (j + 1));

            glm::uvec3 ceCnt{0,0,0};
            glm::uvec3 oe1Cnt{0,0,0};
            unsigned oe2Cntx = 0;
            unsigned oe2Cnty = 0;
            unsigned oe3Cntx = 0;

            std::vector<EvalEdge> EvalEdges = 
            {
                EvalEdge(currentEdge.xIntersects, ceCnt.x,  0, 0,  true),
                EvalEdge(currentEdge.yIntersects, ceCnt.y,  0, 3,  true),
                EvalEdge(currentEdge.zIntersects, ceCnt.z,  0, 8,  true),
                EvalEdge(currentEdge.yIntersects, ceCnt.y,  1, 1,  false),
                EvalEdge(currentEdge.zIntersects, ceCnt.z,  1, 9,  false),
                EvalEdge(otherEdge1.xIntersects,  oe1Cnt.x, 2, 2,  true),
                EvalEdge(otherEdge1.zIntersects,  oe1Cnt.z, 2, 10, true),
                EvalEdge(otherEdge2.xIntersects,  oe2Cntx,  4, 4,  true),
                EvalEdge(otherEdge2.yIntersects,  oe2Cnty,  4, 7,  true),
                EvalEdge(otherEdge1.zIntersects,  oe1Cnt.z, 3, 11, false),
                EvalEdge(otherEdge2.yIntersects,  oe2Cnty,  5, 5,  false),
                EvalEdge(otherEdge3.xIntersects,  oe3Cntx,  6, 6,  true)
            };


            int xL = std::min(currentEdge.xL, std::min(otherEdge1.xL, std::min(otherEdge2.xL, otherEdge3.xL)));
            int xR = std::max(currentEdge.xR, std::max(otherEdge1.xR, std::max(otherEdge2.xR, otherEdge3.xR)));

            int isYEdgeCube = (j == volume.dims.y - 2);
            int isZEdgeCube = (k == volume.dims.z - 2);

            size_t edgeTriId = triangleCount[k * (volume.dims.y - 1) + j];

            for (size_t i = xL; i < xR; i++)
            {
                std::array<size_t, 12> globalIdxs;
                
                int isXEdgeCube = (i == volume.dims.x - 2);

                int gridEdgeCondition = isXEdgeCube | (isYEdgeCube << 1) | (isZEdgeCube << 2);

                int cubeCase = cubeCases[k * (volume.dims.x - 1) * (volume.dims.y - 1) + j * (volume.dims.x - 1) + i];

                if (cubeCase == 0 || cubeCase == 255)
                    continue;

                std::vector<std::pair<glm::vec3, float>> cubeData = getCubeData(volume, i, j, k);

                for (auto ee : EvalEdges)
                {
                    if (config::cubeEdgeCuts[cubeCase][ee.id])
                    {
                        size_t idx = ee.edgeAxisCnt + ee.extraCnt;
                        if ((gridEdgeCondition & ee.condition) == ee.condition)
                            vertices[idx] = interpolate(cubeData, ee.id);
                        globalIdxs[ee.id] = idx;
                        if (ee.increment)
                            ++ee.extraCnt;
                    }
                }

                for (int trid = 0; config::caseTriangles[cubeCase][trid] != -1; trid += 3)
                {
                    triangles[edgeTriId].x = globalIdxs[config::caseTriangles[cubeCase][trid]];
                    triangles[edgeTriId].y = globalIdxs[config::caseTriangles[cubeCase][trid + 1]];
                    triangles[edgeTriId].z = globalIdxs[config::caseTriangles[cubeCase][trid + 2]];
                    ++edgeTriId;
                }
            }
        }
    }
}

uint8_t FlyingEdges::getCubeCase(const std::array<uint8_t, 4>& cubeEdgeCases)
{
    uint8_t cubeCase = 0;
    for (size_t i = 0; i < 4; i++)
    {
        int ccInt0 = (i % 2) == 0 ? (int)std::pow(2, i * 2) : ((int)std::pow(2, i * 2) * 2);
        int ccInt1 = (i % 2) == 0 ? ((int)std::pow(2, i * 2) * 2) : (int)std::pow(2, i * 2);

        if (cubeEdgeCases[i] == 0 || cubeEdgeCases[i] == 2)
            cubeCase |= ccInt0;
        
        if (cubeEdgeCases[i] == 0 || cubeEdgeCases[i] == 1)
            cubeCase |= ccInt1;
    }

    return cubeCase;
}

glm::vec3 FlyingEdges::interpolate(std::vector<std::pair<glm::vec3, float>> cubeData, int edgeIdx)
{
    // Source:
    // https://cs.stackexchange.com/questions/71102/interpolation-on-marching-cubes-algorithm
    std::pair<glm::vec3, float> vertData0 = cubeData.at(config::vertices[edgeIdx][0]);
    std::pair<glm::vec3, float> vertData1 = cubeData.at(config::vertices[edgeIdx][1]);

    if (std::abs(m_isovalue - vertData0.second) < 0.00001)
        return vertData0.first;

    if (std::abs(m_isovalue - vertData1.second) < 0.00001)
        return vertData1.first;

    if (std::abs(vertData0.second - vertData1.second) < 0.00001)
        return vertData0.first;

    float mu = (m_isovalue - vertData0.second) / (vertData1.second - vertData0.second);
    return {
        interp(vertData0.first.x, vertData1.first.x, mu),
        interp(vertData0.first.y, vertData1.first.y, mu),
        interp(vertData0.first.z, vertData1.first.z, mu)
    };
}

float FlyingEdges::interp(float a, float b, float c)
{
    return a + c * (b - a);
}

std::vector<std::pair<glm::vec3, float>> FlyingEdges::getCubeData(const Volume& volume, int i, int j, int k)
{
    std::vector<std::pair<glm::vec3, float>> cubeData;
    for (int z = k; z < k + 2; z++)
    {
        for (int y = j; y < j + 2; y++)
        {
            for (int x = (y == j) ? (i) : (i + 1); y == j ? (x < i + 2) : (x >= i); y == j ? (x++) : (x--))
            {
                glm::vec3 pos{x * m_spacing, y * m_spacing, z * m_spacing};
                float isoVal = (float)volume.at(x, y, z) / 255.f;
                cubeData.push_back({pos, isoVal});
            }
        }
    }
    return cubeData;
}


} // namespace SurfaceExtractor