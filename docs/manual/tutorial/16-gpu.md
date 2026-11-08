# 16. On the GPU

A large grid asks for millions of distances per regrid. The last chapter runs the batch of chapter 5 on the device
with OpenACC and checks it against the host batch.

<<< @/examples/snippets/tut16_gpu-device.f90

- `enter_device` copies the flat arrays the device kernels read (the BVH nodes, the triangle payload, the
  pseudo-normals, the winding-number cache) to the device, once. The surface must be analyzed with the SAH BVH, the
  default: on an octree `enter_device` refuses with a message on stderr and the surface stays on the host. A second
  call is a no-op.
- `distance_many_device` has the arguments and the defaults of `distance_many` (unsigned and squared unless asked).
  Each call moves only the points in and the distances out. Two signs run on the device, `SIGN_PSEUDO_NORMAL` and
  `SIGN_FAST_WINDING_NUMBER`; any other sign algorithm, or a surface not on the device, is refused with a message on
  stderr.
- `exit_device` releases the device copies. `destroy` does not: call it before the surface is destroyed or goes out
  of scope, and before `analyze` or anything else that rebuilds the tree. `move` (chapter 11) refreshes the device
  copy itself.
- `query_many_device`, `box_classify_many_device`, `wall_distance_many_device` and the curve batches follow the same
  pattern.

## Running it

<<< @/examples/output/tut16.ansi{ansi}

This output comes from a gfortran build. Without an OpenACC compiler the `!$acc` directives are comments, so the
"device" batch runs on the host, through the same kernels, and matches the host batch bit for bit: 0 mismatches of
13824 for both signs. The program is still useful there, as a check that the device code path is correct.

![the plane z = -0.026 of the batch, each point's cell shaded by the signed distance distance_many_device returns with the winding-number sign, teal inside and ember outside, the dragon cut by the plane in yellow; none of the 576 values differs from distance_many](../../examples/images/tut16-slice.svg){.fossil-figure}

The figure is drawn from that run: the device code path, executed on the host, not a GPU.

On a GPU, build FOSSIL with nvfortran 26.9 or later (NVIDIA HPC SDK). The `static-nvidia` mode builds the library,
`tests-nvidia` the distance benchmark; the `local_nvf` varset sets the compute capability (`cc89` by default, edit
`$NVF_CC` for another GPU, e.g. `cc80` for an A100):

```bash
fobis build --mode static-nvidia --varset local_nvf
fobis build --mode tests-nvidia --varset local_nvf
nvfortran -acc -gpu=cc89,nofma -Mnofma -I static/mod docs/examples/src/tut16_gpu.f90 static/libfossil.a -o tut16_gpu
```

The flags `-gpu=nofma -Mnofma` turn off fused multiply-add contraction on the device and on the host. FMA contraction
differs between the two, and with it the device results are no longer bit-identical to the host ones.

::: tip What you learned
`enter_device` once, `distance_many_device` per batch, `exit_device` at the end; built with gfortran the same calls
run on the host with identical results. Reference: [GPU distance queries](/guide/advanced/device-distance).
:::

This is the end of the tutorial. The [guide](/guide/) and the [reference](/guide/api/) cover every method and argument.
