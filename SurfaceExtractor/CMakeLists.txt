target_sources(SurfaceExtractor PRIVATE
                                main.cpp
    MarchingCubes.hpp           MarchingCubes.cpp
    MarchingTetrahedra.hpp      MarchingTetrahedra.cpp
    ObjExporter.hpp             ObjExporter.cpp
    Volume.hpp                  
    VolumeLoader.hpp            VolumeLoader.cpp
)
