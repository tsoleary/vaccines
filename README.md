# Vaccines: Vaccination Strategies on Genotype Networks
### Evolutionary Computation Final Project
#### Alex Burnham, & Thomas O'Leary, Blake Williams

#### NOTES MONDAY 12/2/19 
-Blake: run 3 vaccine version on real nets, save to a different filename, gets rid of net size 20/43/52
-Blake: Brute force a couple nets of size 20/43/52/82 w/3 vaccines y-axis=time, x-axis=net size, 1=GA 2=
-Tom: plot mean of ~100 random vaccinations, cut plot off at 1-2yrs, take best 1 of (timeToBestSol) * (popsize) many, repeat for as many GA runs
-Alex: break up figures, fix  real net figure w solutions to be more readable, captions for figs, stats on runs so far, n choose k figure, writing on results, add random solution distributions to plots, ->get 3vacc data when available and incorporate, gets rid of net size 20/43/52


#### NOTES THURSDAY 11/28/19 
Adjusted # functions calls to be up to the point at which the best solution was found.

Generated temporal node/edge files to be used for the growing network part: 'days':= # days since 01/01/2000, indicating for a node when that strain was first sampled, and for an edge the most recent of its two endpoints (ie when that edge comes into existence). This will be used for the growing network to test solutions produced for the older ~1/2 of the network as it grows.

#### NOTES MONDAY 11/25/19 (From the last meeting)
For transcendence=1,2,5, etc
	for each real network of size N
		evolve GA x times 

plot distribution of pct subcritical

Figure 1: Diagram of fitness 

Figure 2: Vaccinations in toy networks (verification/explanatory) 
(a) Star
(b) Chain
(c) Lattice
(d) Erdos Renyi
Show 1 run and explain simple happenings

Table 1:
Parameter

Figure 3: Vaccination coverage vs network size, for transcendence=1,2,3…

Figure 4: Decay of coverage w/network growth (COOL), also w simulated nets

comment on attempts w selection function, mutation function alt

 
#### NOTES SATURDAY 11/23/19
-The framework is there, now its mostly finalizing parameters for 'good' runs, implementing multiple runs (restarts), and making figures.  The Driver's 'finished-ness' decreases as you get further into it.

-Consider playing with the two fitness functions: one determines proportion of nodes 'covered', the other determines mean outbreak size by # nodes (calculates the expected compononent size for an outbreak at a given node). One targets more central strains, the other targets more bridge strains to break up the network.

-Consider trying the '@selectionstochunif' function for selection in place of tournament. I haven't found a reason to use it though.

-Play with transcendence values, slightly, as needed.

-I experimented with a mutation function that would give a ~10% chance for a mutation to pick not randomly from all nodes, but pick from the neighboring nodes of the previous value/node/strain. I thought this would allow solutions to fine-tune themselves as they narrowed in, but I saw minimal change in final fitness/how it got there. If you're curious about this I can put it back in for you to play with. I tested it on the lattice/star/chain and an ER random graph to little effect, but not on the real ones since I didn't see justification for it (but maybe it would offer an improvement there).


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
