/*!
 *  @author    Dubsky Tomas
 */
#pragma once
#include <array>

#include <SurfaceExtractor/Volume.hpp>

namespace SurfaceExtractor {

class MarchingTetrahedra {
public:
    static void extractSurface(
        const Volume& vol,
        std::vector<glm::vec3>& vertices,
        std::vector<glm::uvec3>& triangles
    );
private:
    static void marchTetrahedron(
        const std::array<unsigned int, 4u>& indices,
        const std::array<float, 8u>& corners,
        const glm::vec3& pos,
        std::vector<glm::vec3>& vertices,
        std::vector<glm::uvec3>& triangles
    );
    static const glm::vec3 k_vertexOffsets[8];
    static const int k_indexTable[16][8];
};

}