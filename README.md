Symmetric Molecular Dynamics Engine
=========================

Compiling
-------------------------

### Dependencies

 * libgsl-dev (GNU Scientific Library headers)
 * CMake

```sh
mkdir build && cd build
cmake ..
make && make install
```
The executable is `symd`

Python Interface
-----------------
TODO: write this section

Try running `python Tester.py` in the python directory

Example
-------------------------
Examples may be run with:

    cd example_example
    symd run_params_t.json


Movies
----------

#f5f4e9, (0.9607843137254902, 0.9568627450980393, 0.9137254901960784)


```
~/ffmpeg-git-20210724-amd64-static/ffmpeg -i Trajectories/rods.%05d.bmp -framerate 60 -tune animation -preset slower  -crf 18 -c:v h264 -movflags +faststart -vf "format=yuv420p,setpts=PTS,drawtext=text='@_172135352171_':fontsize=48:x=(w-text_w)/2:y=(h - text_h):fontcolor=#333333:fontfile=/mnt/c/Users/white/Downloads/Courier_Prime/CourierPrime-Regular.ttf" rods.mp4
```
