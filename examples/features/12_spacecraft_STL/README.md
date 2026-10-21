# 12 Spacecraft from STL: 3U CubeSat

**Feature:** complex geometry: STL import (overlapping closed shells handled by generalised winding
numbers) combined with primitives by CSG (thin solar panels, boom), an electrically isolated
probe, mesothermal O+ flow with B and motional E, photoemission with self-shadowing, Sternglass
SEE, per-panel surface currents.

Regenerate the STL with `python3 ../../geometry/make_cubesat_stl.py ../../geometry/cubesat_3u.stl`.
The full run takes about 6 minutes on 2 cores.

## Result

The spacecraft floats at −1.53 V and the isolated probe at −1.53 V. Sunlit faces have a positive
net current density (photoemission exceeds collection), shadowed faces collect electrons, and an
O⁺ wake with downstream focusing forms behind the body.

![render](../../../docs/media/12_cubesat_render.png)

![slices](../../../docs/media/12_cubesat_slices.png)

![wake](../../../docs/media/12_cubesat_wake.gif)
