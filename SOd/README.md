# CEMP for SO(d) synchronization

The metric for CEMP is chosen as the ``Frobenious norm of the difference`` (regard as Stiefel manifold). 

## Usage

Check out the code in ``Examples/Compare_algorithms_SOd.m``, and the following is a sample output under nonuniform and self-consistent corruption.


```    
    Algorithms    MeanError 
    __________    __________

    "Spectral"       0.24842
    "SDP"            0.35469
    "IRLS"         0.0017015
    "CEMP+MST"    1.7033e-08
    "CEMP+GCW"    2.0139e-08
```    
    
