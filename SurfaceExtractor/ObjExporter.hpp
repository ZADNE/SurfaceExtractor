/*!
 *  @author    Dubsky Tomas
 */
#pragma once
#include <string>
#include <vector>

#include <glm/vec3.hpp>

namespace SurfaceExtractor {

class ObjExporter {
public:
    static void save(
        const std::string& outFilename,
        const std::vector<glm::vec3>& vertices,
        const std::vector<glm::uvec3>& triangles
    );
};

}