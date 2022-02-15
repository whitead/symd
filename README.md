Symmetric Molecular Dynamics Engine
=========================

TODO:

Find TODO

Nlist broken?
1. Revise to use enclosing cubes and check implementation
2. Investigate - might just be inconsistent when we had jumps but is ok


NPT
1. volume - https://en.wikipedia.org/wiki/Levi-Civita_symbol#Generalization_to_n_dimensions, https://math.stackexchange.com/questions/1605690/how-to-extend-the-parallelepiped-volume-formula-to-higher-dimensions

3D
1. Write out matrix inverse
2. Compiler Flag

Science
1. Wyckoff Monte Carlo
2. BAOAB
3. 3D


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
~/ffmpeg-git-20210724-amd64-static/ffmpeg -framerate 60 -i Trajectories/cm.%05d.bmp  -crf 18 -c:v h264 -movflags +faststart -vf "format=yuv420p,setpts=PTS,drawtext=text='@_172135352171_':fontsize=48:x=(w-text_w)/2:y=(h - text_h):fontcolor=#333333:fontfile=$HOME/miniconda3/
envs/ternviz/fonts/open-fonts/IBMPlexMono-Light.ttf" traj_cm.mp4
```