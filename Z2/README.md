# Cycle-Edge Message Passing for Z2 Synchronization

## Usage
Checkout and run the following demo code. 
```
Examples/Compare_algorithms_Z2.m
```
The following is a sample output under nonuniform corruption, which shows great advantage of CEMP:

```
Algorithms    MeanError
__________    _________

"Spectral"     0.23483 
"SDP"          0.26782 
"CEMP+GCW"    0.039106 

```

We do not recommend ``CEMP+MST`` for discrete groups and thus only include ``CEMP+GCW``.
