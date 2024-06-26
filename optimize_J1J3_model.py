# pylint: disable=invalid-name
# Colloboration between Franz Bethke and Melanie Koser (Humboldt-Universität zu Berlin)
"""Minimize spin field energy.

    E(u) = -alpha u^T C1 u + u^TC2 u 

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

"""

import os
import sys

from datetime import datetime

from scipy.optimize import minimize  # type: ignore
from scipy.optimize import NonlinearConstraint

from IPython import embed

import matplotlib.pyplot as plt  # type: ignore
import numpy as np
import scipy.sparse as scs  # type: ignore


def get_shift_mat(shape, shift, axis, keep_shape=False):
    """Construct a sparse matrix that represents a shift operation on a tensor.

    Reshaping a tensor to a vector, then applying the shift matrix from the
    left and the reshaping again results in the desired shift of the tensor.

        shifted_u = (D@u.reshape(-1)).reshape(u.shape)
    """

    if shift == 0:
        return scs.eye(np.prod(shape))

    if axis >= len(shape):
        shape += (axis - len(shape) + 1)*[1]

    # the matrix is built by constructing its nontrivial (off)diagonals

    # number of nontrivial diagonals
    ndiags = np.prod(shape[:axis], dtype=int)
    # length of each diagonal
    length = np.prod(shape)
    # 1's in each diagonal
    nentries = np.prod(shape[axis:], dtype=int)

    # number of consecutive trivial rows of the matrix
    nzeros = np.abs(shift)*np.prod(shape[axis+1:], dtype=int)

    diagonals = [np.zeros(length) for _ in range(ndiags)]
    offsets = [-(shift > 0)*nzeros for _ in range(ndiags)]
    keep_indices = []

    for k, diag in enumerate(diagonals):
        offsets[k] -= k*nzeros
        diag[k*nentries:(k+1)*nentries] = 1.
        if shift > 0:
            keep_indices += list(range(k*nzeros + k*nentries,
                                       k*nzeros + (k+1)*nentries))
        else:
            keep_indices += list(range((k+1)*nzeros + k*nentries,
                                       (k+1)*nzeros + (k+1)*nentries))

    out_length = (np.prod(shape[:axis], dtype=int)
                  * (shape[axis] + np.abs(shift))
                  * np.prod(shape[axis+1:], dtype=int))

    mat = scs.diags(diagonals, offsets, (out_length, length), format='csr')

    if keep_shape:
        mat = mat[keep_indices, :]

    return mat


def print_tensor(u):
"""  Function printing the spin field u """
    if len(u.shape) < 3:
        print(u)
        return

    mat = np.empty(u.shape[:-1], dtype=object)

    for i in range(u.shape[0]):
        for j in range(u.shape[1]):
            mat[i, j] = ", ".join([str(u[i, j, k])
                                   for k in range(u.shape[-1])])

    print(mat)


class SpinFieldConfig:
    """ Class of spin fields u, with the attributes 
    delta : interaction parameter, 
    alpha : 4(1-delta) (the original interaction parameter),
    shape : (n, m, 2), where n is the number of rows and m the number of columns, and two the number of entries
            of the vector u_{i,j},
    eps : lattice with, here eps = 1/n, 
    eps_sq : eps^2,
    b : left boundary condition being equal to [0,1],
    Cb : Shift matrix for the boundary layer on the left,
    Cx : Shift matrix for boundary layer on the bottom,
    q : energy on the left boundary layer, 
    H : Internal energy

    Moreover, we introduce 
    C1 : Shift matrix in horizontal direction,
    C2 : Shift matrix in vertical direction, 
    C : matrix for the quadratic spin energy E = U C U^T"""

    def __init__(self, delta, n, m=None):

        self.delta = delta
        self.alpha = 4*(1-delta)

        m = m if m is not None else n-1
        self.shape = (n, m, 2)
        self.eps = 1/n
        self.eps_sq = self.eps**2

        self.b = np.tile(np.array([0, 1.]), (n, 1, 1)).reshape(-1)

        self.Cb = get_shift_mat((n, 1, 2), -m, 1, False)
        self.Cx = get_shift_mat(self.shape, 1, 1, False)

        C1 = get_shift_mat((n, m+1, 2), 1, 1, True)
        C2 = get_shift_mat((n, m+1, 2), 2, 1, True)
        C = -self.alpha*C1 + C2

        R1 = get_shift_mat(self.shape, 1, 0, True)
        R2 = get_shift_mat(self.shape, 2, 0, True)
        R = -self.alpha*R1 + R2

        self.q = self.eps_sq*self.b.T@self.Cb.T@(C + C.T)@self.Cx
        Q = self.eps_sq*(self.Cx.T@C@self.Cx + R)
        self.H = Q + Q.T

    def to_tensor(self, x):
        """ Transforms sparse matrix to tensor """
        return x.reshape(self.shape)

    def get_hspeed(self, x):
        """ Computes the horizontal speed with respect to a spin field 
        Input: self (class of spin configurations)
        Output: Field of angle velocity field """
        bx = self.Cb@self.b + self.Cx@x
        field = bx.reshape(self.shape[0], self.shape[1] + 1, 2)

        C1 = get_shift_mat(field.shape, -1, 1, True)
        C1field = (C1@bx).reshape(field.shape)
        hspeed = np.arccos(np.sum(C1field * field, axis=2))
        hsign = np.where(np.cross(C1field, field) > 0, 1, -1)

        return hsign * hspeed

    def get_vspeed(self, x):
        """ Same as above for the vertical speed """ 
        bx = self.Cb@self.b + self.Cx@x
        field = bx.reshape(self.shape[0], self.shape[1] + 1, 2)

        R1 = get_shift_mat(field.shape, -1, 0, True)
        R1field = (R1@bx).reshape(field.shape)
        vspeed = np.arccos(np.sum(R1field * field, axis=2))
        vsign = np.where(np.cross(R1field, field) > 0, 1, -1)

        return vsign * vspeed

    def plot(self, x, axes=None, cbaxis=None, **kwargs):
        """ Plot the Spin configuration """
        axes = axes if axes is not None else plt.subplots(1, 2,
                                                          figsize=(14, 7))[1]

        bx = self.Cb@self.b + self.Cx@x
        field = bx.reshape(self.shape[0], self.shape[1] + 1, 2)

        hspeed = self.get_hspeed(x)
        vspeed = self.get_vspeed(x)

        curls = (vspeed[:-1, :-1]
                 + hspeed[1:, :-1]
                 - vspeed[:-1, 1:]
                 - hspeed[:-1, :-1])

        cmin = min(np.min(hspeed[:, :-1]), np.min(vspeed[:-1, :]))
        cmax = max(np.max(hspeed[:, :-1]), np.max(vspeed[:-1, :]))
        clim = [max(cmin - np.sign(cmin)*0.2*cmin, -np.pi),
                min(cmax + np.sign(cmax)*0.2*cmax, np.pi)]

        for ax, c, title in zip(axes, [hspeed, vspeed], ['hori', 'vert']):

            ax.set_aspect('equal', 'box')

            if field.shape[0] < 15 and field.shape[1] < 15:
                ax.set_xticks(np.arange(field.shape[1]))
                ax.set_yticks(np.arange(field.shape[0]))
                ax.grid(linestyle=':')

            for i in range(curls.shape[0]):
                for j in range(curls.shape[1]):

                    if -1 < curls[i, j] < 1:
                        continue

                    color = 'red' if curls[i, j] > 1 else 'blue'
                    rect = plt.Rectangle((j, i), 1, 1, color=color, alpha=0.4)
                    ax.add_patch(rect)

            quiver = ax.quiver(field[:, :, 0], field[:, :, 1],
                               c,
                               # cmap='copper',
                               pivot='tail',
                               angles='xy',
                               scale=1.,
                               units='xy',
                               clim=clim,
                               **kwargs)

            ax.set_title(title)

        if cbaxis is None:
            cb = axes[0].figure.colorbar(quiver,
                                         ax=axes.ravel().tolist(),
                                         cax=cbaxis)
                                        # location='top')
            cbaxis = cb.ax
        else:
            cb = plt.colorbar(quiver, cax=cbaxis) #location='top'

         #ticks = cb.get_ticks().tolist()
         #ticks.append(np.arccos(1-self.delta))
         #cb.set_ticks(ticks)
         #ticks[-1] = r'$\theta_{opt}$'
        # cb.set_ticklabels(ticks)

        return axes, cbaxis

    def energy(self, x):
        """ Energy of a Spin configuration """
        return self.q.T@x + 1/2*x.T@self.H@x

    def grad_energy(self, x):
        """ The gradient of the energy """
        return self.q + self.H@x

    def hess_energy(self, _):
         """ The hessian of the energy """
        return self.H

    def constraints(self, x):
         """ Introduce the constraints as an L2 penalization"""
        norms = np.linalg.norm(self.to_tensor(x), axis=2).reshape(-1)
        return 0.5*(norms**2 - 1)

    def jac_constraints(self, x):
        return scs.block_diag(x.reshape(-1, 1, 2), format='csr')


def get_spin_field_guess(spin_config, guess='const rot'):
    """ Input: Spin configuration, guess (choice of start value)
        Output: Spin configuration of the start value depending on the guess:
        
        Ferromagnet : Constant spin field
        const rot: Optimal profile everywhere except on the left boundary
        Laminate : There is one laminate between the different optimal profiles
        Multi-laminate : multiple laminates between the different optimal profiles 
        multi vortex : A vortex construction with vortices close to the left boundary
        random : a random start configuration
        some curls : configuration for testing """

    delta = spin_config.delta
    n, m = spin_config.shape[:-1]

    theta = np.arccos(1-delta)

    # Critical number of vortices
    n_vortex = np.pi // theta

    up = np.array([0, 1.])
    u = np.zeros((n, m, 2))

    for i in range(n):

        for j in range(m):

            if guess == 'ferrormagnet':
                u[i, j, :] = up

            elif guess == 'const rot':

                rot_ij = np.array(
                    [[np.cos((i+j+1) * theta), -np.sin((i+j+1) * theta)],
                     [np.sin((i+j+1) * theta), np.cos((i+j+1) * theta)]])

                # set the value at the point (x, y) = (j*eps, i*eps)
                u[i, j, :] = rot_ij @ up

            elif guess == 'laminate':
                x = j  # *eps
                y = i  # *eps

                if x < y < (n-1)-x:
                    u[i, j, :] = up
                elif x >= y and y <= (n-1)/2:
                    rot_ij = np.array(
                        [[np.cos((-i+j+1) * theta), -np.sin((-i+j+1) * theta)],
                         [np.sin((-i+j+1) * theta), np.cos((-i+j+1) * theta)]])
                    u[i, j, :] = rot_ij @ up
                elif y >= (n-1)-x:
                    rot_ij = np.array(
                        [[np.cos((-(n-i)+j+1) * theta), -np.sin((-(n-i)+j+1) * theta)],
                         [np.sin((-(n-i)+j+1) * theta), np.cos((-(n-i)+j+1) * theta)]])
                    u[i, j, :] = rot_ij @ up
                else:
                    raise AssertionError

            elif guess == 'multi laminate':
                x = j  # *eps, horizontal
                y = i%n_vortex  # *eps, vertical

                if x < y < (n_vortex - 1) - x:
                    u[i, j, :] = up
                elif (x >= y and y <= (n_vortex - 1) / 2) or ((n_vortex -1)/2 < x and y <= (n_vortex - 1) / 2):
                    rot_ij = np.array(
                        [[np.cos((-y + x + 1) * theta), -np.sin((-y + x + 1) * theta)],
                         [np.sin((-y + x + 1) * theta), np.cos((-y + x + 1) * theta)]])
                    u[i, j, :] = rot_ij @ up
                elif y >= (n_vortex - 1) - x or (n_vortex -1)/2 < x:
                    rot_ij = np.array(
                        [[np.cos((-(n_vortex - y) + x + 1) * theta), -np.sin((-(n_vortex - y) + x + 1) * theta)],
                         [np.sin((-(n_vortex - y) + x + 1) * theta), np.cos((-(n_vortex - y) + x + 1) * theta)]])
                    u[i, j, :] = rot_ij @ up
                else:
                    raise AssertionError

            elif guess == 'multi vortex':
                x = j  # *eps, horizontal
                y = i % (2 * n_vortex)  # *eps, vertical

                if x < y  <= n_vortex:
                    u[i, j, :] = up
                elif (x >= y and y <= n_vortex - 1) or n_vortex < x:
                    rot_ij = np.array(
                        [[np.cos((y + x + 1) * theta), -np.sin((-y + x + 1) * theta)],
                         [np.sin((-y + x + 1) * theta), np.cos((-y + x + 1) * theta)]])
                    u[i, j, :] = rot_ij @ up
                elif y >= (n_vortex - 1) and x <= n_vortex: # Construct a vortex
                    y = y % n_vortex
                    if y <= x:
                        rot_ij = np.array(
                            [[np.cos(np.pi + (x-y)*np.pi/x), -np.sin(np.pi + (x-y)*np.pi/x)],
                             [np.sin(np.pi + (x-y)*np.pi/x), np.cos(np.pi + (x-y)*np.pi/x)]])
                        u[i, j, :] = rot_ij @ up
                    else:
                        rot_ij = np.array(
                            [[np.cos(x*np.pi/y), -np.sin(x*np.pi/y)],
                             [np.sin(x*np.pi/y), np.cos(x*np.pi/y)]])
                        u[i, j, :] = rot_ij @ up
                else:
                    raise AssertionError

            elif guess == 'random':
                spin = np.random.randn(2)
                u[i, j, :] = spin / np.linalg.norm(spin)

            elif guess == 'some curls':
                u[i, j, :] = up

                if i == 1 and j == 1:
                    u[i, j, :] = np.array([-1, 0])
                if i == 2 and j == 0:
                    u[i, j, :] = np.array([1, 0])
                if i == 2 and j == 1:
                    u[i, j, :] = np.array([0, -1])

                if i == 5 and j == 0:
                    u[i, j, :] = np.array([1, 0])
                if i == 5 and j == 1:
                    u[i, j, :] = np.array([0, -1])
                if i == 6 and j == 1:
                    u[i, j, :] = np.array([-1, 0])

            else:
                raise ValueError('Unknown guess type!')

    return u.reshape(-1)


def main():

    delta = 0.015
    n = 100
    m = 20

    sfc = SpinFieldConfig(delta, n, m)

    x = get_spin_field_guess(sfc, guess='multi vortex')
    initial_energy = sfc.energy(x)

    dname = datetime.now().strftime('%Y%m%d%H%M') + f'_{delta=}_{n=}_{m=}'
    if not os.path.exists(f'{dname}'):
        os.mkdir(f'{dname}')
    else:
        print('results folder already exists')
        sys.exit(1)

    plt.ion()
    axes, cbaxis = sfc.plot(x)
    fig = axes[0].figure
    axes[0].figure.suptitle(f'it: {0}\nE(u) = {initial_energy}')
    plt.savefig(f'{dname}/inital.svg')

    fig.canvas.draw()
    fig.canvas.flush_events()


    # List of the results
    Steps = [0]
    results=[initial_energy]

    def callback(x, res, axes=axes, cbaxis=cbaxis):
        Steps.append(res.nit)

        np.savetxt(f'{dname}/x_{res.nit}.txt', res.x)

        fig = axes[0].figure

        for ax in fig.axes:
            ax.clear()

        sfc.plot(x, axes, cbaxis)
        axes[0].figure.suptitle(f'it: {res.nit}\nE(u) = {res.fun}')

        fig.canvas.draw()
        fig.canvas.flush_events()
        results.append(sfc.energy(x))

    normed = NonlinearConstraint(sfc.constraints,
                                 lb=0., ub=0.,
                                 jac=sfc.jac_constraints)

    result = minimize(sfc.energy,
                      x,
                      method='trust-constr',
                      jac=sfc.grad_energy,
                      hess=sfc.hess_energy,
                      constraints=[normed],
                      options={'verbose': 2, 'gtol': 1e-6, 'maxiter': 1000},
                      callback=callback,
                      )

    axes, cbaxis = sfc.plot(result.x)
    fig = axes[0].figure
    axes[0].figure.suptitle(f'it: {result.nit}\nE(u) = {sfc.energy(result.x)}')
    plt.savefig(f'{dname}/result.svg')

    print(f"{result = }")
    print(f"{initial_energy = }")

    plt.figure()
    plt.scatter(Steps, results)
    plt.xlabel('Iteration')
   # plt.ylabel('E')
    plt.title('Convergence rate')
    plt.savefig(f'{dname}/convergenc_rate.svg')

    plt.show()
    # embed()



if __name__ == "__main__":
    main()


