/*!
 *  @author    Dubsky Tomas
 */
#include <SurfaceExtractor/MarchingTetrahedra.hpp>

#include <glm/common.hpp>

namespace SurfaceExtractor {

constexpr float k_1over255 = 1.0f / 255.0f;
constexpr float k_threshold = 31.0f / 255.0f;

void MarchingTetrahedra::extractSurface(
    const Volume& vol,
    std::vector<glm::vec3>& vertices,
    std::vector<glm::uvec3>& triangles
) {
    vertices.reserve(4096 * 3);
    triangles.reserve(4096);
    unsigned int zs = vol.dims.x * vol.dims.y;                  //z-stride
    unsigned int ys = vol.dims.y;                               //y-stride

    //For the whole volume
    for (unsigned int z = 0u; z < (vol.dims.z - 1u); ++z) {
        unsigned int zo = zs * z;
        for (unsigned int y = 0u; y < (vol.dims.y - 1u); ++y) {
            unsigned int yzo = zo + ys * y;
            for (unsigned int x = 0u; x < (vol.dims.x - 1u); ++x) {
                glm::vec3 pos{x, y, z};
                unsigned int xyzo = yzo + x;
                std::array<float, 8u> corners;
                unsigned int index = 0u;
                for (unsigned int i = 0u; i < 8u; ++i) {        //Fetch corners of the cube
                    glm::uvec3 offset{
                        (((i + 1) >> 1u) & 1u),
                        ((i >> 1u) & 1u),
                        ((i >> 2u) & 1u)
                    };
                    uint8_t corner =
                        vol.data[xyzo + offset.x + offset.y * ys + offset.z * zs];
                    float val = static_cast<float>(corner) * k_1over255;
                    index |= static_cast<unsigned int>(val > k_threshold) << i;
                    corners[i] = val;
                }
                if (index == 0u || index == 255u) {             //If cube is entirely inside or outside of the volume
                    continue;                                   //Skip this cube because there is no surface to create
                }

                //March the 6 tetrahedra of the cube
                marchTetrahedron({0, 2, 3, 7}, corners, pos, vertices, triangles);
                marchTetrahedron({0, 2, 6, 7}, corners, pos, vertices, triangles);
                marchTetrahedron({0, 4, 6, 7}, corners, pos, vertices, triangles);
                marchTetrahedron({0, 6, 1, 2}, corners, pos, vertices, triangles);
                marchTetrahedron({0, 6, 1, 4}, corners, pos, vertices, triangles);
                marchTetrahedron({5, 6, 1, 4}, corners, pos, vertices, triangles);
            }
        }
    }
}

void MarchingTetrahedra::marchTetrahedron(
    const std::array<unsigned int, 4u>& indices,
    const std::array<float, 8u>& corners,
    const glm::vec3& pos,
    std::vector<glm::vec3>& vertices,
    std::vector<glm::uvec3>& triangles
) {
    unsigned int index = 0u;
    for (unsigned int i = 0u; i < 4u; ++i) {                    //Analyze corners of the tetrahedron
        index |= static_cast<unsigned int>(corners[indices[i]] > k_threshold) << i;
    }
    const int* table = k_indexTable[index];

    auto interpolate = [&](unsigned int i, unsigned int j) {
        float interp = (k_threshold - corners[indices[i]]) / (corners[indices[j]] - corners[indices[i]]);
        if (interp < 0.0f || interp > 1.0) {
            int stop = 5;
        }
        return pos + glm::mix(k_vertexOffsets[indices[i]], k_vertexOffsets[indices[j]], interp);
    };

    if (table[0] != -1) {                                       //Generate first triangle
        unsigned int firstIndex = static_cast<unsigned int>(vertices.size());
        vertices.emplace_back(interpolate(table[0], table[1]));
        auto v1 = interpolate(table[2], table[3]);
        vertices.emplace_back(v1);
        auto v2 = interpolate(table[4], table[5]);
        vertices.emplace_back(v2);
        triangles.emplace_back(firstIndex, firstIndex + 1u, firstIndex + 2u);
        if (table[6] != -1) {                                   //Generate second triangle
            vertices.emplace_back(interpolate(table[6], table[7]));
            triangles.emplace_back(firstIndex + 1u, firstIndex + 2u, firstIndex + 3u);
        }
    }
}

const glm::vec3 MarchingTetrahedra::k_vertexOffsets[8] = {
    {0.0f, 0.0f, 0.0f},
    {1.0f, 0.0f, 0.0f},
    {1.0f, 1.0f, 0.0f},
    {0.0f, 1.0f, 0.0f},
    {0.0f, 0.0f, 1.0f},
    {1.0f, 0.0f, 1.0f},
    {1.0f, 1.0f, 1.0f},
    {0.0f, 1.0f, 1.0f}
};

const int MarchingTetrahedra::k_indexTable[16][8] = {
    {-1, -1, -1, -1, -1, -1, -1, -1},
    { 0,  1,  0,  2,  0,  3, -1, -1},
    { 1,  0,  1,  3,  1,  2, -1, -1},
    { 0,  3,  0,  2,  1,  3,  1,  2},
    { 2,  0,  2,  1,  2,  3, -1, -1},
    { 0,  3,  0,  1,  2,  3,  1,  2},
    { 1,  3,  0,  1,  2,  3,  0,  2},
    { 3,  0,  3,  2,  3,  1, -1, -1},
    { 3,  0,  3,  2,  3,  1, -1, -1},
    { 1,  3,  0,  1,  2,  3,  0,  2},
    { 0,  3,  0,  1,  2,  3,  1,  2},
    { 2,  0,  2,  1,  2,  3, -1, -1},
    { 0,  3,  0,  2,  1,  3,  1,  2},
    { 1,  0,  1,  3,  1,  2, -1, -1},
    { 0,  1,  0,  2,  0,  3, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1}
};

}