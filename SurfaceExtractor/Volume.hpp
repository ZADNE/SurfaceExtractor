/*!
 *  @author    Dubsky Tomas
 */
#pragma once
#include <vector>

#include <glm/vec3.hpp>

namespace SurfaceExtractor {

struct Volume {
    glm::uvec3 dims;
    std::vector<uint8_t> data;

    uint8_t at(unsigned int x, unsigned int y, unsigned int z) const {
        return data[dims.x * dims.y * z + dims.x * y + x];
    }
};

}