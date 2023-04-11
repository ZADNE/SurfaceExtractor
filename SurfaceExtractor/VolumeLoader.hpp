/*!
 *  @author    Dubsky Tomas
 */
#pragma once
#include <string>
#include <vector>

#include <glm/vec3.hpp>

namespace SurfaceExtractor {

struct Volume {
    glm::vec3 dims;
    std::vector<uint8_t> data;
};

class VolumeLoader {
public:
    static Volume load(const std::string& name);
};

}