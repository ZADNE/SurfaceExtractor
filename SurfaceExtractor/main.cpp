/*!
 *  @author    Dubsky Tomas
 */
#include <iostream>

#include <SurfaceExtractor/VolumeLoader.hpp>

int main(int argc, char* argv[]) {
    using namespace SurfaceExtractor;

    try {
        std::string inputFile = (argc >= 2) ? argv[1] : "";

        auto volume = VolumeLoader::load(inputFile);

    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}