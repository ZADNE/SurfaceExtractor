/*!
 *  @author    Dubsky Tomas
 */
#include <SurfaceExtractor/ObjExporter.hpp>

#include <fstream>

namespace SurfaceExtractor {

void ObjExporter::save(
    const std::string& outFilename,
    const std::vector<glm::vec3>& vertices,
    const std::vector<glm::uvec3>& triangles
) {
    std::ofstream f{ outFilename, std::ios::trunc };

    for (const auto& v : vertices) {
        f << "v " << v.x << ' ' << v.y << ' ' << v.z << '\n';
    }

    for (const auto& t : triangles) {
        f << "f " << (t.x + 1) << ' ' << (t.y + 1) << ' ' << (t.z + 1) << '\n';
    }
}

}