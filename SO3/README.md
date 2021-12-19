# CEMP for SO(3) synchronization (Rotation Averaging)

The metric for CEMP is chosen as the geodesic distance on SO(3).

## Usage

Check out the code in ``Examples/Compare_algorithms_SO3.m``, and the following is a sample output under nonuniform and self-consistent corruption.



```
    Algorithms    MeanError 
    __________    __________

    "Spectral"       0.16544
    "SDP"             0.1905
    "IRLS"          0.015304
    "CEMP+MST"    1.2166e-08
    "CEMP+GCW"    0.00089073

```
