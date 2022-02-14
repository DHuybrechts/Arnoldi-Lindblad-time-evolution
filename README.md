# Arnoldi-Lindblad time evolution: Faster-than-the-clock algorithm for the spectrum of time-independent and Floquet open quantum systems

[1] F. Minganti and D. Huybrechts, Arnoldi-Lindblad time evolution: Faster-than-the-clock algorithm for the spectrum of time-independent and Floquet open quantum systems, [Quantum 6, 649](https://doi.org/10.22331/q-2022-02-10-649) (2022).

## Popular summary:

Characterizing the dynamical properties of open quantum systems is a major theoretical challenge. For instance, determining how long quantum properties survive the action of dissipation is pivotal for the development of quantum technologies. Such a task is accomplished by the diagonalization of the Liouvillian superoperator, i.e., the mathematical object encoding both the coherent and dissipative evolution of a system. The Liouvillian eigendecomposition, however, is numerically challenging, and it is common practice to investigate dynamical properties by simply integrating the corresponding equations of motion. These types of extrapolations, however, discard much of the information about the system’s evolution, and become inefficient and imprecise in many cases.

Here, we present a new algorithm that combines the efficiency of the time evolution with the precision and effectiveness of the Liouvillian diagonalization. We perform a series of consecutive short evolutions on appropriately renormalized density matrices, consistent with an Arnoldi iteration algorithm performed on the system’s dynamical map. Accordingly, we construct a reduced Liouvillian incorporating all the dynamical properties otherwise “hidden” by the complexity of the dissipative dynamics. The diagonalization of this much smaller reduced Liouvillian retrieves the long-living dynamical properties up to numerical precision and, moreover, manages to determine the long-time dynamics from a much shorter one, making it a faster-than-the-clock algorithm. Our algorithm is general and can be applied to any system described by a Lindblad master equation, including Floquet (i.e., periodically driven) ones.

## Article abstract:

The characterization of open quantum systems is a central and recurring problem for the development of quantum technologies. For time-independent systems, an (often unique) steady state describes the average physics once all the transient processes have faded out, but interesting quantum properties can emerge at intermediate timescales. Given a Lindblad master equation, these properties are encoded in the spectrum of the Liouvillian whose diagonalization, however, is a challenge even for small-size quantum systems. Here, we propose a new method to efficiently provide the Liouvillian spectral decomposition. We call this method an Arnoldi-Lindblad time evolution, because it exploits the algebraic properties of the Liouvillian superoperator to efficiently construct a basis for the Arnoldi iteration problem.
The advantage of our method is double: (i) It provides a faster-than-the-clock method to efficiently obtain the steady state, meaning that it produces the steady state through time evolution shorter than needed for the system to reach stationarity. (ii) It retrieves the low-lying spectral properties of the Liouvillian with a minimal overhead, allowing to determine both which quantum properties emerge and for how long they can be observed in a system. This method is *general and model-independent*, and lends itself to the study of large systems where the determination of the Liouvillian spectrum can be numerically demanding but the time evolution of the density matrix is still doable. Our results can be extended to time evolution with a time-dependent Liouvillian. In particular, our method works for Floquet (i.e., periodically driven) systems, where it allows not only to construct the Floquet map for the slow-decaying processes, but also to retrieve the stroboscopic steady state and the eigenspectrum of the Floquet map. Although the method can be applied to any Lindbladian evolution (spin, fermions, bosons, ...), for the sake of simplicity we demonstrate the efficiency of our method on several examples of coupled bosonic resonators (as a particular example). Our method outperforms other diagonalization techniques and retrieves the Liouvillian low-lying spectrum even for system sizes for which it would be impossible to perform exact diagonalization.


## The example code:

The function for the Arnoldi-Lindblad time evolution [1] can be found in `ArnoldiLindbladTimeEvolution.py` and has been applied to an example in the Jupyter notebook file `ArnoldiLindblad_example.ipynb`.
The function uses QuTiP [2, 3] for time evolution of the density matrix. If one wishes to use the example function without QuTiP one has to change the “mesolve” function to your own time evolution function.


## Future outlook:

An application of the method will be published on the arXiv in the foreseeable future (and updated here) where we utilise the Arnoldi-Lindblad time evolution to calculate the low-lying Liouvillian spectrum in a permutational invariant system consisting of N n-level systems [4]. It allows to calculate the eigenvalues of a Liouvillian whose matrix representation in the full space of states is equal to n^(2N) x n^(2N), with actual dimensions e.g. of the order 7.10^(10) x 7.10^(10) , far beyond the capabilities of exact diagonalization. We thus combine the advantages of a permutational invariant procedure with the Arnoldi-Lindblad time evolution and will verify the (non-)occurence of time crystalline behaviour that is predicted by mean-field theory.


## References:

[1] F. Minganti and D. Huybrechts, Arnoldi-Lindblad time evolution: Faster-than-the-clock algorithm for the spectrum of time-independent and Floquet open quantum systems, [Quantum 6, 649](https://doi.org/10.22331/q-2022-02-10-649) (2022).

[2] J. Johansson, P. Nation and F. Nori, QuTiP: An open-source Python framework for the dynamics of open quantum systems, [Comput. Phys. Commun. 183, 1760](https://www.sciencedirect.com/science/article/abs/pii/S0010465512000835?via%3Dihub) (2012).

[3] J. Johansson, P. Nation and F. Nori, QuTiP 2: A Python framework for the dynamics of open quantum systems, [Comput. Phys. Commun. 184, 1234](https://www.sciencedirect.com/science/article/abs/pii/S0010465512003955?via%3Dihub) (2013)

[4] To be published (in preparation)
