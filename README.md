## Surface Extractor
Implementation of three methods for extracting surfaces from volumetric data. This project is implemented as part of the VGE course on VUT FIT, Brno.

Authors:
Boris Burkalo,
Tomáš Dubský

This solution contains implementation for the Marching cubes, Marching tetrahedra and Flying edges algorithms.

## Downloading the data
This sample code directly works with specified data which can be downloaded [here](https://drive.google.com/drive/folders/1ZodXP7_f4KIrUSYIN80-SaBAr2F6fxjn?usp=share_link). Please download the .zip file and extract it within the root of this repository.

## Running the program
When the data is successfully placed into the root directory just run:
```
mkdir build && cd build
cmake ..
make
```
The application then can be found in the `build` directory and can be ran as:
```
./SurfaceExt -i [filepath] -o [filename].obj -a [ c | e | t ] -iso [(0.0, 1.0)] 
```
where:
 - `-i [filepath]` is a filename to one of the test data within build repository
 - `-o [filename].obj` sets the output `.obj` file
 - `-a [ c | e | t ]` sets the algorithm to be used
   - `c` for Marching cubes
   - `e` for Flying edges
   - `t` for Marching tetrahedra
 - `-iso [(0.0, 1.0)]` is a isosurface threshold within the interval 0.0 and 1.0