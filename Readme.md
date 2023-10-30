### Based on :
Eigen  
Spectra  
MMG  

### File contents in TestOn2D:
| FileName | Contents |
|-----:|-----------|
|E_worst_K_inverse.cpp| The maximum eigenvalue of the inverse of the stiffness matrix K is directly calculated.|
|E_worst_EP_Pp.cpp| The eigenvalue problem(EP), which is the method 3 mentioned in the previous PDF. (Error) |
|E_worst_GEP_Pp.cpp|The broad eigenvalue problem is computed using eigen, where P is defined as $P=HH^T$, which is the basis vector matrix of p |