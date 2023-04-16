/*!
 *  @author    Dubsky Tomas
 */
#pragma once
#include <SurfaceExtractor/Volume.hpp>

namespace SurfaceExtractor {

class MarchingCubes {
public:
    static void extractSurface(
        const Volume& vol,
        std::vector<glm::vec3>& vertices,
        std::vector<glm::uvec3>& triangles
    );
private:
    static const int k_adjacentVertexIndices[12];
    static const glm::vec3 k_vertexOffsets[8];
    static const int k_indexTable[256][16];
};

}