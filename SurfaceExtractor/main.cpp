/*!
 *  @author    Dubsky Tomas
 */
#include <iostream>

#include <SurfaceExtractor/VolumeLoader.hpp>
#include <SurfaceExtractor/ObjExporter.hpp>

namespace SurfaceExtractor {

enum class ExtractionAlgorithm {
    MARCHING_CUBES,
    FLYING_EDGES
};

bool equals(const char* str, std::initializer_list<const char*> options) {
    for (const auto& opt : options) {
        if (std::strcmp(str, opt) == 0) {
            return true;
        }
    }
    return false;
}

}

int main(int argc, char* argv[]) {
    using namespace SurfaceExtractor;

    try {
        //Parse command line arguments
        std::string inputFile = "teapot.raw";
        std::string outputFile = "out.obj";
        ExtractionAlgorithm algorithm = ExtractionAlgorithm::MARCHING_CUBES;
        for (int i = 1; i < (argc - 1); i += 2) {
            if (equals(argv[i], {"-a", "--algorithm"})) {
                if (equals(argv[i + 1], {"c", "cubes", "marching_cubes"})) {
                    algorithm = ExtractionAlgorithm::MARCHING_CUBES;
                } else if (equals(argv[i + 1], {"e", "edges", "flying_edges"})) {
                    algorithm = ExtractionAlgorithm::FLYING_EDGES;
                } else {
                    std::cerr << "Unknown algorithm \"" << argv[i + 1] << "\"\n";
                }
            } else if (equals(argv[i], {"-i", "--input"})) {
                inputFile = argv[i + 1];
            } else if (equals(argv[i], {"-o", "--output"})) {
                outputFile = argv[i + 1];
            } else {
                std::cerr << "Unknown option \"" << argv[i + 1] << "\"\n";
            }
        }

        //Load input
        auto volume = VolumeLoader::load(inputFile);

        //Extract the surface
        std::vector<glm::vec3> vertices;
        std::vector<glm::uvec3> triangles;
        switch (algorithm) {
        case ExtractionAlgorithm::MARCHING_CUBES:
            /* Magic goes here... */
            break;
        case ExtractionAlgorithm::FLYING_EDGES:
            /* Magic goes here... */
            break;
        default: throw std::exception{"No algorithm selected"};
        }

        //Save the output
        ObjExporter::save(outputFile, vertices, triangles);

    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}