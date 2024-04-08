**Optimization of the J1-J3 model** 

Our goal is to optimize a special discrete lattice interaction energy, the J1-J3 model on a planar rectangle, e.g., $\Omega \coloneqq [0,1)^2$. 
Let $\alpha \in (0, 4)$ be an interaction parameter describing the strength of the ferromagnetic and antiferromagnetic interaction.
This parameter decides how strong the ferromagnetic and  anti-ferromagnetic interactions are. We assume the atoms of our
material are ordered within a regular square lattice with the lattice width $\varepsilon \in (0,1/2)$. 
To each lattice point $(i\varepsilon, j\varepsilon) \in \Omega\cap\varepsilon\mathbb{Z}^2$ we associate a normed vector $u_{i,j} \in \mathbb{S}^1$, the so-called Spin vector.
We want to understand the properties of the ground states of the energy

   $$ I(u) \coloneqq -\alpha \sum_{i,j} u_{i,j} \cdot u_{i+1, j} - \alpha \sum_{i,j} u_{i,j} \cdot, u_{i, j+1} + \sum_{i,j} u_{i,j} \cdot u_{i+2, j} + \sum_{i,j} \cdot u_{i, j} \cdot u_{i, j+2}$$
    
for a spin field $u \colon \Omega \cap \varepsilon \mathbb{Z}^2\rightarrow \mathbb{S}^1$ with $u_{0, \cdot} = (0,1)^T$.

**For an overview work, we refer to:**

  Marco Cicalese, Marwin Forster, and Gianluca Orlando. “Variational analysis of a two-
  dimensional frustrated spin system: emergence and rigidity of chirality transitions”.

  Hung T. Diep, ed. Frustrated spin systems. 2nd ed.


**Packages needed:**

  os,
  sys,
  datetime,
  scipy,
  IPython,
  matplotlib,
  numpy

**Optimization with Scipy from a chosen start value**:

*Ferromagnet* : Constant spin field

*Const rot*: Optimal profile everywhere except on the left boundary

*Laminate* : There is one laminate between the different optimal profiles

*Multi-laminate* : Multiple laminates between the different optimal profiles 

*Multi vortex* : A vortex construction with vortices close to the left boundary

*Random* : A random start configuration

*Some curls* : Vortex configuration for testing


**Output:** Generates figures of the optimized spin field during the optimization. The intermediate
 optimized spin fields are saved in an associated folder. Furthermore, this folder will contain
 svg images of the start configuration and the end configuration.

**Expected energy regimes:**
        
Let $\delta \coloneqq (4-\alpha)/4\in (0,1)$ the given material parameter in the following. We expect 
the following regimes:

  1. **Ferromagnet** if $\delta^{1/2} \leq \varepsilon$:
     Constant spin field with energy $\delta^2$ and no vortices

  2.  **Helimagnet** if $\delta^{1/2} \exp(-1/\delta) < \varepsilon < \delta^{1/2}$:
      Multi-laminate more or less with energy
      $\varepsilon\delta^{3/2} (\vert ln(\varepsilon/\delta^{1/2})\vert +1)$ with no vortices

  3. **Vortex structure** if $\varepsilon < \delta^{1/2} \exp(-1/\delta)$: Spin field using
      vortices with the energy $\varepsilon \delta^{1/2}$

For more details we refer to a upcoming work of Ginster, Koser, and Zwicknagl 2024.

*Credits to Franz Bethke and Melanie Koser*



