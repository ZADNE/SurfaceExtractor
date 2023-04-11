/*!
 *  @author    Dubsky Tomas
 */
#include <iostream>

#include <SurfaceExtractor/VolumeLoader.hpp>
#include <SurfaceExtractor/ObjExporter.hpp>

int main(int argc, char* argv[]) {
    using namespace SurfaceExtractor;

    try {
        //Load input
        std::string inputFile = (argc >= 2) ? argv[1] : "teapot.raw";
        auto volume = VolumeLoader::load(inputFile);

        //Extract the surface
        std::vector<glm::vec3> vertices;
        std::vector<glm::uvec3> triangles;
        /* Magic goes here... */

        //Save the output
        std::string outputFile = (argc >= 3) ? argv[2] : "out.obj";
        ObjExporter::save(outputFile, vertices, triangles);

    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}