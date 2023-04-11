/*!
 *  @author    Dubsky Tomas
 */
#include <SurfaceExtractor/VolumeLoader.hpp>

#include <fstream>
#include <array>

namespace SurfaceExtractor {

constexpr auto VOLUME_NAMES = std::to_array<const char*>({
    "foot.raw",
    "bonsai.raw",
    "teapot.raw",
    "head.raw",
    "manix.raw",
    "tooth.raw"
});
constexpr auto VOLUME_DIMS = std::to_array<glm::uvec3>({
    {256u, 256u, 256u},
    {256u, 256u, 256u},
    {256u, 256u, 178u},
    {256u, 256u, 225u},
    {128u, 128u, 115u},
    {256u, 256u, 151u}
});
static_assert(VOLUME_NAMES.size() == VOLUME_DIMS.size());


Volume VolumeLoader::load(const std::string& filename) {
    //Select file to load
    std::string volName;
    Volume vol;
    for (size_t i = 0; i < VOLUME_NAMES.size(); i++) {
        if (filename == VOLUME_NAMES[i]) {
            volName = VOLUME_NAMES[i];
            vol.dims = VOLUME_DIMS[i];
            break;
        }
    }
    if (volName == "") {
        throw std::exception{"Unknown input filename"};
    }

    //Load the volume
    std::ifstream file{volName, std::ios::binary | std::ios::ate};
    std::streamsize size = file.tellg();
    if (size < (vol.dims.x * vol.dims.y * vol.dims.z)) {
        throw std::exception{"Input file is smaller than expected"};
    }
    file.seekg(0, std::ios::beg);
    vol.data.resize(size);
    static_assert(sizeof(char) == sizeof(uint8_t));
    file.read(reinterpret_cast<char*>(vol.data.data()), size);

    return vol;
}

}