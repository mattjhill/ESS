ESS
===

Extended SemiSeparable solver

To run the solver key in the following:

	make -f makefile.mk

Then run the executable with the number of unknowns and rank of the sem-separability

	./ESS N m c

where `N` is the number of unknowns, `m` is the rank of semi-separability and `c` is a character to denote, if we want to compare with usual method or not. For instance,

	./ESS 5000 4 y

means you are a solving a `5000` by `5000` dense system with semi-separable rank being `4` and would like to compare with the usual dense solver.

	./ESS 500000 4 n

means you are a solving a `500000` by `500000` dense system with semi-separable rank being `4` and would not like to compare with usual dense solver. Typically, larger the value of `N` (say beyond `10000`), it is not advisable to compare with the usual solver, since the usual solver will take excruciatingly long time to solve.