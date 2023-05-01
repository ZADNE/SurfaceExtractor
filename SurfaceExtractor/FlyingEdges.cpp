#include <SurfaceExtractor/FlyingEdges.hpp>
#include <SurfaceExtractor/Volume.hpp>
#include <SurfaceExtractor/FlyingEdgesCases.hpp>

#include <iostream>
#include <cmath>
#include <limits>
#include "FlyingEdges.hpp"

#define SPACING 1.f

namespace SurfaceExtractor
{

float FlyingEdges::m_isovalue = 31.f / 255.f;

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
                    if (j != volume.dims.y - 1)
                    {
                        size_t curIdY = k * volume.dims.y * (volume.dims.x - 1) + (j + 1) * (volume.dims.x - 1) + i;
                        uint8_t curEdgeCaseY = edgeCases[curIdY];
                        
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

            for (size_t i = xL; i < xR; i++)
            {
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

                // Edges that belong to the row.
                currentEdge.xIntersects += config::cubeEdgeCuts[cubeCase][0];
                currentEdge.yIntersects += config::cubeEdgeCuts[cubeCase][3];
                currentEdge.zIntersects += config::cubeEdgeCuts[cubeCase][8];

                // Add the values from cubes that are on the edge of the
                // voxel grid.
                bool xEnd = i == (volume.dims.x - 2);
                bool yEnd = j == (volume.dims.y - 2);
                bool zEnd = k == (volume.dims.z - 2);
                if (xEnd)
                {
                    currentEdge.yIntersects += config::cubeEdgeCuts[cubeCase][1];
                    currentEdge.zIntersects += config::cubeEdgeCuts[cubeCase][9];
                }

                if (yEnd)
                {
                    otherEdge1.xIntersects += config::cubeEdgeCuts[cubeCase][2];
                    otherEdge1.zIntersects += config::cubeEdgeCuts[cubeCase][10];
                }

                if (zEnd)
                {
                    otherEdge2.xIntersects += config::cubeEdgeCuts[cubeCase][4];
                    otherEdge2.yIntersects += config::cubeEdgeCuts[cubeCase][7];
                }

                if (xEnd && yEnd)
                    otherEdge1.zIntersects += config::cubeEdgeCuts[cubeCase][11];
                
                if (xEnd && zEnd)
                    otherEdge2.yIntersects += config::cubeEdgeCuts[cubeCase][5];

                if (yEnd && zEnd)
                    otherEdge3.xIntersects += config::cubeEdgeCuts[cubeCase][6];
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
            glm::uvec3 ceCnt{0,0,0};

            const Edge& otherEdge1 = edges.at(k * volume.dims.y + (j + 1));
            const Edge& otherEdge2 = edges.at((k + 1) * volume.dims.y + j);
            const Edge& otherEdge3 = edges.at((k + 1) * volume.dims.y + (j + 1));

            glm::uvec3 oe1Cnt{0,0,0};
            glm::uvec3 oe2Cnt{0,0,0};
            unsigned oe3Cnt = 0;

            int xL = std::min(currentEdge.xL, std::min(otherEdge1.xL, std::min(otherEdge2.xL, otherEdge3.xL)));
            int xR = std::max(currentEdge.xR, std::max(otherEdge1.xR, std::max(otherEdge2.xR, otherEdge3.xR)));

            bool isYEdgeCube = (j == volume.dims.y - 2);
            bool isZEdgeCube = (k == volume.dims.z - 2);

            size_t edgeTriId = triangleCount[k * (volume.dims.y - 1) + j];

            for (size_t i = xL; i < xR; i++)
            {
                std::array<size_t, 12> globalIdxs;
                
                bool isXEdgeCube = (i == volume.dims.x - 2);

                int cubeCase = cubeCases[k * (volume.dims.x - 1) * (volume.dims.y - 1) + j * (volume.dims.x - 1) + i];

                if (cubeCase == 0 || cubeCase == 255)
                    continue;

                std::vector<std::pair<glm::vec3, float>> cubeData = getCubeData(volume, i, j, k);
                
                if (config::cubeEdgeCuts[cubeCase][0])
                {
                    size_t idx = currentEdge.xIntersects + ceCnt.x;
                    vertices[idx] = interpolate(cubeData, 0);
                    globalIdxs[0] = idx;
                    ++ceCnt.x;
                }

                if (config::cubeEdgeCuts[cubeCase][3])
                {
                    size_t idx = currentEdge.yIntersects + ceCnt.y;
                    vertices[idx] = interpolate(cubeData, 3);
                    globalIdxs[3] = idx;
                    ++ceCnt.y;
                }

                if (config::cubeEdgeCuts[cubeCase][8])
                {
                    size_t idx = currentEdge.zIntersects + ceCnt.z;
                    vertices[idx] = interpolate(cubeData, 8);
                    globalIdxs[8] = idx;
                    ++ceCnt.z;
                }

                // edge cases

                if (config::cubeEdgeCuts[cubeCase][1])
                {
                    if (edgeTriId == 20)
                        std::cout << "CECNT " << ceCnt.y << std::endl;
                    size_t idx = currentEdge.yIntersects + ceCnt.y;
                    if (isXEdgeCube)
                    {
                        vertices[idx] = interpolate(cubeData, 1);
                    }
                    globalIdxs[1] = idx;
                }

                if (config::cubeEdgeCuts[cubeCase][9])
                {
                    size_t idx = currentEdge.zIntersects + ceCnt.z;
                    if (isXEdgeCube)
                    {
                        vertices[idx] = interpolate(cubeData, 9);
                    }
                    globalIdxs[9] = idx;
                }

                if (config::cubeEdgeCuts[cubeCase][2])
                {
                    size_t idx = otherEdge1.xIntersects + oe1Cnt.x;
                    if (isYEdgeCube)
                    {
                        vertices[idx] = interpolate(cubeData, 2);
                    }
                    globalIdxs[2] = idx;
                    ++oe1Cnt.x;
                }

                if (config::cubeEdgeCuts[cubeCase][10])
                {
                    size_t idx = otherEdge1.zIntersects + oe1Cnt.z;
                    if (isYEdgeCube)
                    {
                        vertices[idx] = interpolate(cubeData, 10);
                    }
                    globalIdxs[10] = idx;
                    ++oe1Cnt.z;
                }

                if (config::cubeEdgeCuts[cubeCase][4])
                {
                    size_t idx = otherEdge2.xIntersects + oe2Cnt.x;
                    if (isZEdgeCube)
                    {
                        vertices[idx] = interpolate(cubeData, 4);
                    }
                    globalIdxs[4] = idx;
                    ++oe2Cnt.x;
                }

                if (config::cubeEdgeCuts[cubeCase][7])
                {
                    size_t idx = otherEdge2.yIntersects + oe2Cnt.y;
                    if (isZEdgeCube)
                    {
                        vertices[idx] = interpolate(cubeData, 7);
                    }
                    globalIdxs[7] = idx;
                    ++oe2Cnt.y;
                }

                if (config::cubeEdgeCuts[cubeCase][11])
                {
                    size_t idx = otherEdge1.zIntersects + oe1Cnt.z;
                    if (isXEdgeCube && isYEdgeCube)
                    {
                        vertices[idx] = interpolate(cubeData, 11);
                    }
                    globalIdxs[11] = idx;
                }

                if (config::cubeEdgeCuts[cubeCase][5])
                {
                    size_t idx = otherEdge2.yIntersects + oe2Cnt.y;
                    if (isXEdgeCube && isZEdgeCube)
                    {
                        vertices[idx] = interpolate(cubeData, 5);
                    }
                    globalIdxs[5] = idx;
                }

                if (config::cubeEdgeCuts[cubeCase][6])
                {
                    size_t idx = otherEdge3.xIntersects + oe3Cnt;
                    if (isYEdgeCube && isZEdgeCube)
                    {
                        vertices[idx] = interpolate(cubeData, 6);
                    }
                    globalIdxs[6] = idx;
                    ++oe3Cnt;
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
                glm::vec3 pos{x * SPACING, y * SPACING, z * SPACING};
                float isoVal = (float)volume.at(x, y, z) / 255.f;
                cubeData.push_back({pos, isoVal});
            }
        }
    }
    return cubeData;
}


} // namespace SurfaceExtractor