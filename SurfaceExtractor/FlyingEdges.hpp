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


#pragma once
#include <SurfaceExtractor/Volume.hpp>
#include <array>
#include <utility>

namespace SurfaceExtractor {

struct Edge
{
    int xL;
    int xR;

    unsigned xIntersects;
    unsigned yIntersects;
    unsigned zIntersects;
};

struct EvalEdge
{
    EvalEdge(unsigned eAC, unsigned& ec, int cond, int idd, bool ic) 
    : edgeAxisCnt(eAC), extraCnt(ec), condition(cond),
    id(idd), increment(ic)
    {}

    unsigned edgeAxisCnt;
    unsigned& extraCnt;
    int condition;
    int id;
    bool increment;
};

class FlyingEdges
{

public:
    static void extractSurface(
        const Volume& vol,
        std::vector<glm::vec3>& vertices,
        std::vector<glm::uvec3>& triangles,
        float isoValue = 0.121f
    );

    static float m_isovalue;
    static float m_spacing;

private:
    static void pass1(const Volume& volume, std::vector<uint8_t>& edgeCases,
        std::vector<Edge>& edges);
    static void pass2(const Volume& volume, std::vector<uint8_t>& edgeCases,
        std::vector<Edge>& edges, std::vector<uint8_t>& cubeCases, std::vector<unsigned>& triangleCount);
    static void pass3(const Volume& volume, std::vector<Edge>& edges, std::vector<unsigned>& triangleCount,
        unsigned& trianglesTot, unsigned& verticesTot);
    static void pass4(const Volume& volume, const std::vector<Edge>& edges, const std::vector<unsigned>& triangleCount,
        const std::vector<uint8_t>& cubeCases, std::vector<glm::vec3>& vertices, std::vector<glm::uvec3>& triangles);

    static uint8_t getCubeCase(const std::array<uint8_t, 4>& cubeEdgeCases);
    static glm::vec3 interpolate(std::vector<std::pair<glm::vec3, float>> cubeData, int edgeIdx);
    static float interp(float a, float b, float c);
    static std::vector<std::pair<glm::vec3, float>> getCubeData(const Volume& volume, int i, int j, int k);
    
};

} // namespace SurfaceExtractor