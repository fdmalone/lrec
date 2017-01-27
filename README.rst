LREC
====

Perform Liovillian Recursion to solve for the time evolution of a particular
quantum-mechanical operator. The recursion can be carried out using different inner
products, either the infinte temperature version or the ground state found from exact
diagonalisation (obviously restricted to small system sizes). The resulting basis vectors
of operators can also be output to a file so that the inner product can be calculated
using a higher level method, such as quantum Monte Carlo.

Only works for the one-dimensional Heisenberg model, for lattices of some multiple of 8 in
length (I know madness).

There is very little documentation, and parts of the code are a complete mess. I am mostly
dumping this here in the very slim chance an appendix of my thesis ever gets revisited.
