/*!
 *  @author    Dubsky Tomas
 */
#include <iostream>
#include <cstring>
#include <stdexcept>

#include <SurfaceExtractor/VolumeLoader.hpp>
#include <SurfaceExtractor/ObjExporter.hpp>
#include <SurfaceExtractor/MarchingCubes.hpp>
#include <SurfaceExtractor/MarchingTetrahedra.hpp>
#include <SurfaceExtractor/FlyingEdges.hpp>

namespace SurfaceExtractor {

enum class ExtractionAlgorithm {
    MarchingCubes,
    MarchingTetrahedra,
    FlyingEdges
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
        // Parse command line arguments
        std::string inputFile = "teapot.raw";
        std::string outputFile = "out.obj";
        float isovalue = 0.121f;
        ExtractionAlgorithm algorithm = ExtractionAlgorithm::MarchingCubes;
        for (int i = 1; i < (argc - 1); i += 2) {
            if (equals(argv[i], {"-a", "--algorithm"})) {
                if (equals(argv[i + 1], {"c", "cubes", "marching_cubes"})) {
                    algorithm = ExtractionAlgorithm::MarchingCubes;
                } else if (equals(argv[i + 1], {"t", "tetrahedra", "marching_tetrahedra"})) {
                    algorithm = ExtractionAlgorithm::MarchingTetrahedra;
                } else if (equals(argv[i + 1], {"e", "edges", "flying_edges"})) {
                    algorithm = ExtractionAlgorithm::FlyingEdges;
                } else {
                    std::cerr << "Unknown algorithm \"" << argv[i + 1] << "\"\n";
                }
            } else if (equals(argv[i], {"-iso", "--isovalue"}))
            {
                size_t sz;
                isovalue = std::stof(argv[i + 1], &sz);
                if (isovalue < 0 || isovalue > 1.f)
                {
                    std::cout << "Isovalue must be within [0.0f, 1.0f]." << std::endl;
                    return 1;
                }
            } else if (equals(argv[i], {"-i", "--input"})) {
                inputFile = argv[i + 1];
            } else if (equals(argv[i], {"-o", "--output"})) {
                outputFile = argv[i + 1];
            } else {
                std::cerr << "Unknown option \"" << argv[i + 1] << "\"\n";
            }
        }

        // Load input
        auto volume = VolumeLoader::load(inputFile);

        // Extract the surface
        std::vector<glm::vec3> vertices;
        std::vector<glm::uvec3> triangles;
        switch (algorithm) {
        case ExtractionAlgorithm::MarchingCubes:
            MarchingCubes::extractSurface(volume, vertices, triangles, isovalue);
            break;
        case ExtractionAlgorithm::MarchingTetrahedra:
            MarchingTetrahedra::extractSurface(volume, vertices, triangles, isovalue);
            break;
        case ExtractionAlgorithm::FlyingEdges:
            FlyingEdges::extractSurface(volume, vertices, triangles, isovalue);
            break;
        default: throw std::runtime_error("No algorithm selected");
        }

        std::cout << "Dataset: "                    << inputFile        << std::endl;
        std::cout << "Isosurface value threshold: " << isovalue         << std::endl;
        std::cout << "Generated triangles: "        << triangles.size() << std::endl;
        std::cout << "Generated vertices: "         << vertices.size()  << std::endl;
        std::cout << "Output file: "                << outputFile       << std::endl;

        // Save the output
        ObjExporter::save(outputFile, vertices, triangles);

    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}