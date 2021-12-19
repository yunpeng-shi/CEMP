# Cycle-Edge Message Passing for Angular (SO(2)) Synchronization

The metric for CEMP is chosen as the geodesic distance on U(1).

## Usage

Check out the code in Examples/Compare_algorithms.m, and the following is the sample output.


```
    Algorithms    MeanError
    __________    _________

    "Spectral"     0.20936 
    "SDP"          0.26136 
    "IRLS"         0.17556 
    "CEMP+MST"    0.086138 
    "CEMP+GCW"    0.088543 

```
