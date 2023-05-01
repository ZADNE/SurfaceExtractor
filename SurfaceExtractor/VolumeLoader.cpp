/*!
 *  @author    Dubsky Tomas
 */
#include <SurfaceExtractor/VolumeLoader.hpp>

#include <fstream>
#include <array>
#include <stdexcept>

namespace SurfaceExtractor {

constexpr auto k_volumeNames = std::to_array<const char*>({
    "foot.raw",
    "bonsai.raw",
    "teapot.raw",
    "head.raw",
    "manix.raw",
    "tooth.raw"
});
constexpr auto k_volumeDims = std::to_array<glm::uvec3>({
    {256u, 256u, 256u},
    {256u, 256u, 256u},
    {256u, 256u, 178u},
    {256u, 256u, 225u},
    {128u, 128u, 115u},
    {256u, 256u, 151u}
});
static_assert(k_volumeNames.size() == k_volumeDims.size());


Volume VolumeLoader::load(const std::string& filename) {
    //Select file to load
    std::string volName;
    Volume vol;
    for (size_t i = 0; i < k_volumeNames.size(); i++) {
        if (filename == k_volumeNames[i]) {
            volName = k_volumeNames[i];
            vol.dims = k_volumeDims[i];
            break;
        }
    }
    if (volName == "") {
        throw std::runtime_error{"Unknown input filename"};
    }

    //Load the volume
    std::ifstream file{volName, std::ios::binary | std::ios::ate};
    std::streamsize size = file.tellg();
    if (size < (vol.dims.x * vol.dims.y * vol.dims.z)) {
        throw std::runtime_error{"Input file is smaller than expected"};
    }
    file.seekg(0, std::ios::beg);
    vol.data.resize(size);
    static_assert(sizeof(char) == sizeof(uint8_t));
    file.read(reinterpret_cast<char*>(vol.data.data()), size);

    return vol;
}

}