# knot-engine
Generates a Wavefront .OBJ file for a 3D Lissajous knot.
Just run `python3 ./knot-engine.py` and see what happens, or
`python3 ./knot-engine.py --help` for more options.

I find it nice to view them in [MeshLab](https://www.meshlab.net/).
Here's the default output:

![default-curve](screenshots/knot-default.png)

You can make it look janky by setting a low `--curve-density`.
This one is done with `--curve-density 30`:

![janky-curve](screenshots/knot-curve-density-30.png)

The parametric formulae that specify the curve contain three coefficients, given in
`-a`, `-b`, and `c`. This curve was generated with ``-a 1 -b 2 -c 3``:

![curve-1-2-3](screenshots/knot-1-2-3.png)

Here's the same curve with ``--wrap-radius 50``:

![curve-50](screenshots/knot-1-2-3-50.png)
