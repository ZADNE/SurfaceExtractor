/*!
 *  @author    Dubsky Tomas
 */
#pragma once
#include <string>

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
    static const int INDEX_TABLE[256][16];
};

}