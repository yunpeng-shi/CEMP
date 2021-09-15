# Cycle-Edge Message Passing for Z2 Synchronization

## Usage
Download matlab files to the same directory. Checkout and run the following demo code. 
```
Examples/Compare_algorithms_Z2.m
```
The following is a sample output (error in degrees) under Nonuniform and Adversarial corruption, which shows great advantage of CEMP and MPLS:

```
Algorithms    MeanError
__________    _________

"Spectral"     0.23483 
"SDP"          0.26782 
"CEMP+GCW"    0.039106 

```

Here ``Spectral`` and ``SDP`` refers to eigenvector method and SDP relaxation method for approximately solving least squares formulation. 
See [Angular Synchronization by Eigenvectors and Semidefinite Programming,](https://arxiv.org/abs/0905.3174) Amit Singer, Applied and Computational Harmonic Analysis, 2011 for details.
We do not recommend ``CEMP+MST`` for discrete groups and thus only include ``CEMP+GCW``.
