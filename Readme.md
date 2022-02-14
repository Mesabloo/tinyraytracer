# Homework assignment

Final work~:

- Normal rendering:
![](https://raw.githubusercontent.com/mesabloo/tinyraytracer/main/out_normal.jpg)

- Parallax rendering:
![](https://raw.githubusercontent.com/mesabloo/tinyraytracer/main/out_parallax.jpg)

- Stereoscope rendering:
![](https://raw.githubusercontent.com/mesabloo/tinyraytracer/main/out_stereo.jpg)

## compilation
```sh
git clone --recurse-submodules https://github.com/mesabloo/tinyraytracer.git
cd tinyraytracer
git checkout homework_assignment
git submodule update --init
mkdir build
cd build
cmake ..  
make
```
