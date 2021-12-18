# Robust Group Synchronization via Cycle-Edge Message Passing (CEMP)

## Introduction and References

Cycle-edge message passing (CEMP) is a very useful algorithm for robust group synchronization (GS). Examples of GS problem include ``correlation clustering`` (Z2 group), ``phase synchronization`` (SO(2) group), ``rotation averaging`` (SO(3) group), and ``multi-object matching`` (Sn group).

The GS problem asks to recover group elements <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{g_{i}}"> from their noisy/corrupted relative measurements <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{g_{ij}"> (ideally equals to <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{g_{i}g_{j}^{-1}}">).

[Robust Group Synchronization via Cycle-Edge Message Passing](https://link.springer.com/content/pdf/10.1007/s10208-021-09532-w.pdf), Gilad Lerman and Yunpeng Shi, Foundations of Computational Mathematics, 2021.

If you would like to use our code for your paper, please cite the above work.



Here ``Spectral`` refers to eigenvector method for approximately solving least squares formulation. See [Angular Synchronization by Eigenvectors and Semidefinite Programming,](https://arxiv.org/abs/0905.3174) Amit Singer, Applied and Computational Harmonic Analysis, 2011 for details.

``IRLS-GM`` refers to iteratively reweighted least squares (IRLS) that uses Geman-McClure loss function. See [Efficient and Robust Large-Scale Rotation Averaging](https://www.cv-foundation.org/openaccess/content_iccv_2013/papers/Chatterjee_Efficient_and_Robust_2013_ICCV_paper.pdf), Avishek Chatterjee and Venu Madhav Govindu, ICCV 2013 for details.

``IRLS-L0.5`` refers to iteratively reweighted least squares (IRLS) that uses L_1/2 loss function. See [Robust Relative Rotation Averaging](http://www.ee.iisc.ac.in/labs/cvl/papers/robustrelrotavg.pdf), Avishek Chatterjee and Venu Madhav Govindu, TPAMI, 2018 for details.

``CEMP+MST``, ``CEMP+GCW``, ``MPLS`` refer to our methods. ``CEMP+GCW`` only appears in CEMP paper [Robust Group Synchronization via Cycle-Edge Message Passing](https://link.springer.com/content/pdf/10.1007/s10208-021-09532-w.pdf), Gilad Lerman and Yunpeng Shi, Foundations of Computational Mathematics, 2021. It combines CEMP with a weighted spectral method (using graph connection weight matrix).

## A Variety of Corruption Models
We provide 5 different corruption models. 3 for nonuniform topology and 2 for uniform toplogy (see ``Uniform_Topology.m`` and ``Nonuniform_Topology.m``). Uniform/Nonuniform toplogy refers to whether the corrupted subgraph is Erdos Renyi or not. In other words, the choice of Uniform/Nonuniform toplogy decides how to select edges for corruption. In ``Uniform_Topology.m``, two nodes are connected with probability ``p``. Then edges are independently drawn with probability ``q`` for corruption. In ``Nonuniform_Topology.m``, two nodes are connected with probability ``p``. Then with probability ``p_node_crpt`` a node is selected so that its neighboring edges are candidates for corruption. Next, for each selected node, with probability ``p_edge_crpt`` an edge (among the neighboring edges of the selected node) is corrupted. This is a more malicious scenario where corrupted edges have cluster behavior (so local coruption level can be extremely high). 

One can also optionally add noise to inlier graph and outlier graph through setting ``sigma_in`` and ``sigma_out`` for ``Nonuniform_Topology.m``. For ``Uniform_Topology.m`` we assume inlier and outlier subgraph have the same noise level ``sigma``.

The argument ``crpt_type`` in the two functions determines how the corrupted relative rotations are generated for those selected edges. In ``Uniform_Topology.m``, there are 2 options of ``crpt_type``: ``uniform`` and ``self-consistent``.
In ``Nonuniform_Topology.m``, there are the following 3 options of ``crpt_type``.

``uniform``: The corrupted relative permutations <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{g_{ij}}"> are i.i.d follows uniform distribution over the space of SO(3).

``self-consistent``: The corrupted <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{g_{ij}}"> are relative rotations of another set of absolute rotations. Namely <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{g_{ij} = g_i^{crpt} g_j^{crpt}'}"> where those absolute rotations are different from the ground truth and are i.i.d drawn from the uniform distribution in the space of SO(3). In this way, the corrupted relative rotations are also cycle-consistent.

``adv``: Extremely malicious corruption that replaces the underlying absolute rotations from ground truth <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{g_i^*}"> to <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{g_i^{crpt}}">. Namely <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{g_{ij} = g_i^{crpt} g_j^{* }'}"> for the corrupted neighboring edges (i,j) of node i. Additional high noise must be added to the outlier-subgraph, otherwise the recovery of the ground truth can be ill-posed. This model is not included in MPLS paper. It was first introdued in [Robust Multi-object Matching via Iterative Reweighting of the Graph Connection Laplacian, NeurIPS 2020](https://proceedings.neurips.cc/paper/2020/file/ae06fbdc519bddaa88aa1b24bace4500-Paper.pdf) for permutation synchronization.



## Implementation of CEMP

The demo code uses the following function for implementing CEMP:
```
CEMP(Ind, RijMat, parameters)
```
Each row of ``Ind`` matrix is an edge index (i,j). The edge indices (the rows of Ind) MUST be sorted in ``row-major order``. That is, the edge indices are sorted as  for example (1,2), (1,3), (1,4),..., (2,3), (2,5), (2,8),..., otherwise the code may crash when some edges are not contained in any 3-cycles. Make sure that i<j. If some edges have indices (3,1), then change it to (1,3) and take a transpose to the corresponding Rij.


