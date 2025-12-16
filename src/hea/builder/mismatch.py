import math
import numpy as np
from inspect import currentframe
from ase import Atoms

def cfrac_approx(x,kmax=None,nmin=None,nmax=None,include_next=False,verbose=False):
    '''
    x =approx= h/k

    modified from: https://www.embeddedrelated.com/showarticle/1620.php

    input:
    kmax : max value of k
    nmin : min number of iterations (has priority over kmax, to avoid 1/1)
    nmax : max number of iterations
    include_next : ?

    returns:
    # coeffs : integers of continued fraction (dont really care)
    brow[1] : h
    brow[2] : k
    '''
    if kmax is None and nmax is None:
        raise ValueError('Need to specify either max denominator kmax or max number of coefficients nmax')
    if kmax is None:
        kmax = float('inf')
    if nmax is None:
        nmax = float('inf')
    if nmin is None:
        nmin=3
    arow = [1,0,1]
    brow = [x,1,0]
    coeffs = []
    n = 0
    while n < nmax:
        q = math.floor(x)
        r = x - q
        if r == 0:
            break
        n += 1
        x = 1/r
        xrow = [x, arow[1]+q*brow[1], arow[2]+q*brow[2]]
        if n>nmin:
           #if (xrow[2] > kmax and not include_next) or brow[2] > kmax:
           if (xrow[2] > kmax and not include_next):
              break
        arow = brow
        brow = xrow
        coeffs.append(q)
        if verbose:
            print(brow, q)
    if len(coeffs) > 1 and coeffs[-1] == 1:
        coeffs = coeffs[:-1]
        coeffs[-1] += 1
    # return coeffs, brow[1], brow[2]
    return brow[1], brow[2]



def lat_mismatch( a0_a, brav_a, n_layers_a, \
                  a0_b, brav_b, n_layers_b, \
                  interlattice_dist=None, vacuum=None, \
                  nmin=None, kmax=None, verbose=True ):
    """
    Construct a periodic cell with two lattice types A, and B, where A is on bottom,
    B is on top, and we control the amount of vacuum and the interlattice dist.
    The supercell size is computed such to minimize the mismatch:
    a0_a/a0_b ~= h/k where h and k are integers.
    Works only for cubic brav lattices atm.
    The final conf should look like:

          ______________
          |            |  \
          |            |  |--> :: vacuum
          |            |  /
          |  B  B  B  B|  \
          | B  B  B  B |  |--> :: brav_b * n_layers_b
          |  B  B  B  B|  |
          | B  B  B  B |  /
          |            |  --> :: interlattice_dist
          |A A A A A A |  \
          | A A A A A A|  |--> :: brav_a * n_layers_a
          |A A A A A A |  /
          |____________|
           <- k*a0_a -> :: supercell dimension on x, and y

    **== input: ==**

    :param a0_a: lattice constant of system A
    :type a0_a: float

    :param brav_a: bravais lattice type of system A (possible: 'bcc')
    :type brav_a: str

    :param n_layers_a: number of unit cell repetitions of A in z-direction
    :type n_layers_a: int

    :param a0_b: lattice constant of system B
    :type a0_b: float

    :param brav_b: bravais lattice type of system B (possible: 'bcc')
    :type brav_b: str

    :param n_layers_b: number of unit cell repetitions of B in z-direction
    :type n_layers_b: int

    :param interlattice_dist: distance between z-values of topmost atms of A, and lowermost atoms of B,
                              defulat==0.0 (atoms of B start at last value of A, i.e. they overlap).
                              This is in units of distance (angstrom), projected on z ax.
    :type interlattice_dist: float

    :param vacuum: how much vacuum to put on top of B (empty space in z-direction), units of distance (Ang.)
    :type vacuum: float

    :param nmin: minimal number of iterations for the fractional algo (default=3), this has priority over
                 `kmax`, to avoid cases of 1/1 when `kmax` is too small.
    :type nmin: int

    :param kmax: maximal value of k for the fraction
    :type kmax: int

    :param verbose: write some additional details about boxsize (default=True)
    :type verbose: logical

    **== Output: ==**

    :param atm: object with integer types (`numbers`), positions, and cell.
    :type atm: ase.Atoms object
    """
    if interlattice_dist is None:
        interlattice_dist = 0.0
        print( "Warn at:",currentframe().f_code.co_filename,"line:",currentframe().f_lineno,":"\
              " `interlattice_dist` is zero, atoms might be overlapping.")
    if vacuum is None:
        vacuum = 0.0
        print( "Warn at:",currentframe().f_code.co_filename,"line:",currentframe().f_lineno,":"\
              " `vacuum` is zero, check periodicity in z direction.")


        # default for the fractional algo
    if nmin is None:
        nmin = 3
    if kmax is None:
        kmax= 50

    ## sys A bravais lattice
    if( brav_a == "bcc" ):
        a1_a = np.array([1.0, 0.0, 0.0])*a0_a
        a2_a = np.array([0.0, 1.0, 0.0])*a0_a
        a3_a = np.array([0.0, 0.0, 1.0])*a0_a
        lat_a = np.array([a1_a, a2_a, a3_a]).T
        nat_a = 2
        bas_a = np.array([ \
                           [ 0.0, 0.0, 0.0 ], \
                           [ 0.5, 0.5, 0.5 ] \
                          ])
    else:
        raise ValueError("unsupported value for brav_a:",brav_a)
    # atom types of A
    typ_a = np.ndarray([nat_a], dtype=object)
    typ_a[:] = 1

    ## sys B bravais lattice
    if( brav_b == "bcc" ):
        a1_b = np.array([1.0, 0.0, 0.0])*a0_b
        a2_b = np.array([0.0, 1.0, 0.0])*a0_b
        a3_b = np.array([0.0, 0.0, 1.0])*a0_b
        lat_b = np.array([a1_b, a2_b, a3_b]).T
        nat_b = 2
        bas_b = np.array([ \
                           [ 0.0, 0.0, 0.0 ], \
                           [ 0.5, 0.5, 0.5 ] \
                          ])
    else:
        raise ValueError("unsupported value for brav_b:",brav_b)
    # atom types of B
    typ_b = np.ndarray([nat_b], dtype=object)
    typ_b[:] = 1

    # compute x = h/k
    h, k = cfrac_approx( a0_a/a0_b, nmin=nmin, kmax=kmax )

    if( verbose ):
        print( "frac:",a0_a/a0_b, "has approc integer fractional:",h,"/",k)
        print( "k*a0_a=", k, a0_a, k*a0_a )
        print( "h*a0_b=", h, a0_b, h*a0_b )
        print( "absolute error:", abs(k*a0_a - h*a0_b))
        print( "relative error:", abs(k*a0_a - h*a0_b)/(k*a0_a))

    # shift for sys B, such that it starts at same value of z as atom of A with largest z
    maxz = np.max(bas_b[:,2])
    dz = np.array([0.0, 0.0, n_layers_a-maxz])
    dz = np.matmul( lat_a, dz )

    # total number of atoms
    nat_tot = nat_a*k*k*n_layers_a + nat_b*h*h*n_layers_b
    # total box vectors
    lat_tot = np.array([a1_a*k, \
                        a2_a*k, \
                        a3_a*n_layers_a + a3_b*n_layers_b - \
                        np.matmul(lat_a, np.array([0.0,0.0,maxz])) + \
                        np.array([0.0, 0.0, interlattice_dist+vacuum]) \
                        ])

    # accumulate positions and atomic types
    pos = np.ndarray([nat_tot,3],dtype=float)
    typ = np.ndarray([nat_tot], dtype=int)
    idx=0
    # print( nat_tot )
    # print( 'Lattice="',*lat_tot[0], *lat_tot[1], *lat_tot[2],'"')
    for xx in range(k):
        for yy in range(k):
            for zz in range( n_layers_a):
                shift = np.array([xx,yy,zz],dtype=float)
                for i in range(nat_a):
                    pos[idx] = np.matmul( lat_a, bas_a[i]+shift)
                    typ[idx] = typ_a[i]
                    idx+=1
    for xx in range(h):
        for yy in range(h):
            for zz in range( n_layers_b ):
                shift = np.array([xx,yy,zz],dtype=float)
                for i in range(nat_b):
                    pos[idx] = np.matmul( lat_b, bas_b[i]+shift) + dz+np.array([0.0,0.0,interlattice_dist])
                    typ[idx] = typ_b[i]
                    idx+=1
    # for i in range(nat_tot):
    #     print( typ[i], *pos[i] )
    atm = Atoms(numbers=typ, positions=pos, cell=lat_tot, pbc=True)
    return atm



# test main func
if __name__ == "__main__":
    struc = lat_mismatch( a0_a=3.12, brav_a="bcc", n_layers_a=6, \
                          a0_b=3.19, brav_b="bcc", n_layers_b=9, \
                          interlattice_dist=2.0, vacuum=10.0 \
                   )
