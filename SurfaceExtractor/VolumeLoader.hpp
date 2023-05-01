/*!
 *  @author    Dubsky Tomas
 */
#pragma once
#include <string>

#include <SurfaceExtractor/Volume.hpp>

namespace SurfaceExtractor {

class VolumeLoader {
public:
    static Volume load(const std::string& filename);
};

}