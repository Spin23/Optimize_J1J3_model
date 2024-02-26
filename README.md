Optimization of the J1-J3 model 

Our goal is to optimize a special Heisenberg model, the J1-J3 energy on a planar rectangle, e.g., (0,1)^2. 
Let alpha be an interaction parameter in (0, 4). This parameter decides how strong the ferromagnetic and 
anti-ferromagnetic interactions are. We assume the atoms of our material are ordered within a regular square 
lattice with the lattice width eps contained in (0,1). To each lattice point we associate a normed vector, the 
so-called Spin vector. Let n be the largest integer smaller than 1/eps. We want to understand the properties 
of the ground states of the energy
    E(u) = -alpha u^T C1 u + u^TC2 u,
where u is a n times n dimensional tensor mapping each element to a Spin vector of a lattice point.  

For an overview work we refer to:

  Marco Cicalese, Marwin Forster, and Gianluca Orlando. “Variational analysis of a two-
  dimensional frustrated spin system: emergence and rigidity of chirality transitions”.

  Hung T. Diep, ed. Frustrated spin systems. 2nd ed.

Conventions for variable names:

  eps : the distance between horizontally or vertically neighbouring atoms

  n : the number of rows in the atom grid
        indexed by i

  m : the number of columns in the atom grid
        indexed by j

  u : a matrix of shape n x m x 2
        represents a spin field
        u[i, j] is the 2d vector at (x, y) = (j*eps, i*eps)

  alpha : interaction parameter

  C1 : Shift matrix in horizontal direction

  C2 : Shift matrix in vertical direction

  x : a flattend version of u
        x == u.reshape(-1) should always be true

  Note: A shift of the indices in the first axis (rows) of u by +1, i.e.,
          u[i, j] -> u[i+1, j], ends up beeing an upwards vertical shift of the
          represented spin field y -> y + eps.

   A shift of the indices in the second axis (columns) of u by +1, i.e.,
   u[i, j] -> u[i, j+1], ends up beeing a right horizontal shift of the
   represented spin field x -> x + eps.

  Expected energy regimes:
        In the follwowing let eps = 1/n and delta the given material parameter.

  1. Ferromagnet if delta**0.5 < eps:
     Constant spin field with energy delta**2 and no vortices

  2.  Helimagnet if delta**0.5 * exp(-1/delta) < eps < delta**0.5:
      (Multi-laminate more or less) with energy
      eps*delta^{3/2} (|ln(eps/delta^{1/2})| +1) with no vortices

  3. Vortexmagnet if eps < delta**0.5 * exp(-1/delta):
     energy = eps*delta^{1/2}
     with opt/(eps*2pi) vortices where opt= arccos(1-delta), which is
     the optimal angle velocity.
