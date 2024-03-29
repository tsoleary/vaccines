# Vaccines: Vaccination Strategies on Genotype Networks
### Evolutionary Computation Final Project
#### Blake Williams, Alex Burnham, & Thomas O'Leary

#### Original Proposal
The goal of this project is to determine strain-specific vaccination strategies for a multi-strain disease with an underlying genotype network. A specified number of strains (e.g. 3) will be selected for vaccination, from an underlying genotype network of strains, in a well-mixed population of persons to be infected. The notion is that vaccinating against certain nodes in a genotype network will exploit the benefits of strain-transcending immunity, reducing overall infections by protecting against similar strains. The goal is to evolve and identify the configurations of good vaccination implementations (good= low mean steady state infection counts), as well as the features of nodes found in the best solutions. 

The genome will be represented with a bit-string vector of length n for the n nodes (strains)
of the genotype network. Values of 1 indicate that a strain will be vaccinated against, while a 0
indicates nobody in the well-mixed population is vaccinated against that strain. A simple GA will be
implemented based on the bit-string representation of the genome. Fitness will be evaluated by the
steady-state infection count of an SIRS model run on each vaccination configuration: lower infection
counts will indicate better vaccine distribution configurations. A penalty could be introduced for
increasing the number of vaccine strains in a given genome beyond some threshold. Community
detection may be implemented to improve GA performance with ordering of sites under the suspicion
it will improve crossover, since vaccination strains are expected to be ineffective if they are too close
together.

Code that runs SIRS simulations on genotype networks already exists, which will have to
be modified to include vaccination prior to incorporation with a GA. We will have to decide on test
networks (e.g. lattices and chains) and real-world networks for testing, as well as crossover and
mutation restrictions. We will observe scalability of results with network size, due to superlinear
complexity of the fitness function (SIRS model) that will restrict network size to roughly n < 100,
which should be sufficient to explore the effects of real-world complexity.
