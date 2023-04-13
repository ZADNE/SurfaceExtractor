/*!
 *  @author    Dubsky Tomas
 */
#pragma once
#include <string>
#include <array>

#include <SurfaceExtractor/Volume.hpp>

namespace SurfaceExtractor {

class MarchingCubesExtractor {
public:
    static void extractSurface(
        const Volume& vol,
        std::vector<glm::vec3>& vertices,
        std::vector<glm::uvec3>& triangles
    );
private:
    static void generateTriangles(
        const std::array<float, 8>& corners,
        unsigned int index,
        std::vector<glm::vec3>& vertices,
        std::vector<glm::uvec3>& triangles
    );

    static const int ADJACENT_VERTEX_INDICES[12];
    static const glm::vec3 VERTEX_OFFSETS[8];
    static const int INDEX_TABLE[256][16];
};

}