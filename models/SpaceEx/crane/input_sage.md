```python
sage: from polyhedron_tools.misc import BoxInfty, polyhedron_to_Hrep
sage: B = BoxInfty(lengths=[[0, 5],[0, 0],[-0.2, 0.2],[-0.1, 0.1], [0, 0], [0, 0]])
sage: F, g = polyhedron_to_Hrep(B)
sage: savemat("X0.mat", {"F":F, "g":g})
 ```
 
 ```python
sage: B1 = BoxInfty(lengths=[[0, 5], [0, 0]])
sage: B2 = BoxInfty(lengths=[[-0.2, 0.2], [-0.1, 0.1]])
sage: B3 = BoxInfty(lengths=[[0.0, 0.0], [0.0, 0.0]])
sage: F1, g1 = polyhedron_to_Hrep(B1)
sage: F2, g2 = polyhedron_to_Hrep(B2)
sage: F3, g3 = polyhedron_to_Hrep(B3)
```

```python
sage: from scipy.io import savemat
sage: savemat("X0.mat", {"F":F, "g":g, "F1":F1, "F2":F2, "F3":F3, "g1":g1, "g2":g2, "g3":g3)
```