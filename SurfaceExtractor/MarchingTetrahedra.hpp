/*!
 *  @author    Dubsky Tomas
 */
#pragma once
#include <SurfaceExtractor/Volume.hpp>

namespace SurfaceExtractor {

class MarchingTetrahedra {
public:
    static void extractSurface(
        const Volume& vol,
        std::vector<glm::vec3>& vertices,
        std::vector<glm::uvec3>& triangles
    );
};

}