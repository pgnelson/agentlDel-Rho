import os
import sys
import random 
import math
import time
import pickle
##########################################################################################
#	Name: recordStatistics
#
#	Assumptions: This assumes that the population list's entries will each have 5
#		     sub entries. The 5 sub entries are as following: the rho value of 
#		     the individual as a float, the lDel gene of the individual as a list of
#		     ones and zeroes, the fitness of the individual, the alpha gene of the 
#		     individual, and the beta gene of the individual. It is assumed that they
#	             appear in the order of the last sentences
#
#	Purpose: This function exists for the purpose finding the mean fitness, number of 
#	 	 ones in an lDel gene, and mean rho of an individual in a population. 
#		 Furthermore, the variance of each one of the aforementioned means is 
#		 calculated. Once these statistics are calculated, they are separated by
#		 tabs and written as a line and are written to the outFile object. 
#
#	Arguments: - outFile is file object representing the file where statistics will be 
#		     written to.
#		   - population is a list representing the population that the statistics
#		     will be about
#
#	Returns: Nothing
#
##########################################################################################
def recordStatistics(outFile, population):
	# The statistics to be recorded are initialized to 0
	populationSize = len(population)
	meanFitness = 0
	meanFitnessVariance = 0
	meanRho = 0
	meanlDelLoci = 0
	meanRhoVariance = 0
	meanlDelLociVariance = 0
	meanAlpha = 0
	meanAlphaVariance = 0
	meanBeta = 0
	meanBetaVariance = 0
	
	# Iterate through the population once to find the mean
	for i in range (populationSize):
		meanRho += population[i][0] / float(populationSize)
		meanlDelLoci += population[i][1].count(1) / float(populationSize)
		meanFitness += population[i][2] / float(populationSize)
		meanAlpha += sum(population[i][3]) / float(populationSize)
		meanBeta += sum(population[i][4]) / float(populationSize)
		
	# Iterate another time to find the variance using the means
	for i in range (populationSize):
		meanRhoVariance += (population[i][0] - meanRho)**2 / float(populationSize)
		meanlDelLociVariance += (population[i][1].count(1) - meanlDelLoci)**2 / \
		                        float(populationSize)
		meanFitnessVariance += (population[i][2] - meanFitness)**2 / float(populationSize)
		meanAlphaVariance += (sum(population[i][3]) - meanAlpha)**2 / \
		                          float(populationSize)
		meanBetaVariance += (sum(population[i][4]) - meanBeta)**2 / \
		                          float(populationSize)
	
	# Write the numbers into a text using the join method. 
	numbers = [meanFitness, meanFitnessVariance, meanlDelLoci, meanlDelLociVariance, 
	           meanRho, meanRhoVariance, meanAlpha, meanAlphaVariance, meanBeta, 
	           meanBetaVariance]
	data = [str(statistic) for statistic in numbers]
	outFile.write("\t".join(data) + "\n")
	outFile.flush()
	
##########################################################################################
#	Name: findRhoMaximizingFitness
#
#	Assumptions: When calculating the fitness given a rho, the environmental optimum is 
#		     assumed to be zero. The reason for this is that this function is called
#		     at the population is being initialized and the environmental optimum 
#		     is 0 and the alpha and beta genes' sums are also 0. 
#
#	Purpose: This function exists so that, given an lDel gene, the user can find the rho
#	         value that maximizes the fitness of an individual. This is done by starting
#		 rho at .2. Rho is then multiplied by .99 until it is less than 10**-13, For 
#		 each one of these intermediate rho values, the fitness is calculated with the
#		 the given lDel gene. The highest fitness and the rho corresponding to it are
#		 stored.
#
#	Arguments: s - an integer corresponding to the fitness cost associated with 
#		       proofreading
#	           lDelGene - a list of 1's and 0's whose optimal rho value will be found
#		   pDel -  the probability of a loci in an lDel gene going from benign to 
#			   deleterious. This is only used to calculate the fitness.
#		   pMinusDel - the probability of a loci in lDel gene going from deleterious
#			       to benign. This is only used to calculate fitness.
#			       mulDel - the probability of a given loci in an lDel gene mutating. This is
#						passed to calculate fitness.
#		   alphaGene - the alpha gene that will be used to calculate the fitness
#	           betaGene - the beta gene used to calculate the fitness 
#
#	Returns: A float value in the interval (10**-13, 2] is returned. This float 
#		 corresponds to the rho value that would maximize fitness given the lDel gene
#		 passed by the user. 
#
##########################################################################################
def findRhoMaximizingFitness(s, lDelGene, pDel, pMinusDel, mulDel, alphaGene, betaGene):
	# Start r1ho at .2
	currentMaxFitness = 0
	rhoCorrespondingToMaxFit = 0
	lDelGeneLength = len(lDelGene)
	rho = .2
	
	# Multiply rho by .99  until rho < 10**-13
	while rho > 10**-13:
		# Find the fitness given the current rho and the assumption that the environmental
		# optimum is 0
		fitness = getFitness(s, lDelGene, pDel, pMinusDel, mulDel, rho, alphaGene, 
		                     betaGene, 0)
		                     
	 	# Keep track of the highest fitness and the rho corresponding to it
		if fitness > currentMaxFitness:
			currentMaxFitness = fitness
			rhoCorrespondingToMaxFit = rho
		
		# Increment rho
		rho *= .99
		
	# Return the rho corresponding to the maximum fitness found
	return(rhoCorrespondingToMaxFit)
	
##########################################################################################
#	Name: getFitness
#
#	Assumptions: None
#
#	Purpose: getFitness exists so that one can find the fitness given a collection of 
#		 genetic information, the most important being an lDel gene, a rho value,
#		 an alpha gene and a beta gene.  There are four components of fitness. The
#		 first three components are calculated according the equations in the 2011
#		 Rajon-Masel paper. The last is calculated from the alpha and beta gene beta
#		 gene. A sum over the two lists will be calculated an used as an
#		 environmental readiness that will be compared to the optimum. The sum 
#		 includes alpha genes and possibly beta genes; whether or not these beta
#		 genes are included depends on a strip of lDels at the end of the lDel gene. 
#			
#	Arguments: s - an integer corresponding to the fitness cost associated with 
#		       proofreading
#		   lDelGene - a list of 1's and 0's whose optimal rho value will be found
#		   pDel -  the probability of a loci in an lDel gene going from benign to 
#			   deleterious.
#		   pMinusDel - the probability of a loci in lDel gene going from deleterious
#			       to benign.
#		   mulDel - the probability of a given loci in an lDel gene mutating.
#		   rho - the read through rate of a stop codon 
#		   alphaGene - a list containing ten floats representing 10 different 
#			       pre-stop codon trait values of an individuals
#		   betaGene - a list containing ten floats representing 10 different
#			      post-stop codon trait values
#		   envOptimum - the optimum environmental trait value for an individual
#
#	Returns: A float corresponding to the fitness of an individual is returned. It will
#		 will be in the interval [0, 1] and is the product of the four components of 
#		 fitness.
#
##########################################################################################
def getFitness(s, lDelGene, pDel, pMinusDel, mulDel, rho, alphaGene, betaGene, envOptimum):
	# Calculate the first three the fitness components
	lDelGeneLength = len(lDelGene)
	tempDelFitness = max(0, 1 - s * ((rho *  lDelGene.count(1) / float(lDelGeneLength)) + \
	                     (1 - lDelGene.count(1) / float(lDelGeneLength)) * rho**2 * \
	                     pDel / (pDel + pMinusDel)))		
	permDelFitness = (1 - 23/9 * mulDel)**lDelGene.count(1)
	proofFitness =  1 / (1 - math.log(rho) * 10**-2.5)
	
	# Calculate the environmental fitness
	envTrait = 0
	for i in range(alphaGeneLength):
		expressionOflDel = (lDelGene[len(lDelGene) - (alphaGeneLength - i)] == 1)
		if expressionOflDel:
			envTrait += alphaGene[i]
		else:
			envTrait += alphaGene[i] + rho * betaGene[i]
	envFitness = math.exp(-((envTrait - envOptimum)**2 / (2 * .5**2)))
	
	# Multiply all the components together and return the result 
	fitness = tempDelFitness * permDelFitness * proofFitness * envFitness
	return fitness

##########################################################################################
#	Name: initializePopulation
#
#	Assumptions: It is assumed that all of the individual in the population being 
#		     initialized should have the rho value that makes their fitness the 
#		     highest. It is also assumed the hat each individual has five pieces 
#		     of information the that define them: a rho gene, an lDel, a fitness, a 
#		     beta gene, and an alpha gene.
#
#	Purpose: This function exists to assign genes and a fitness value to each individual 
#		 in a population of a desired size. The genes of an individual include:
#                rho, lDel, alpha, beta. lDel is a list of 1's and 0's, rho is a float, and 
# 		 alpha and beta are lists of float values. An initial proportion for 1's
# 		 is given by the user, and each individual gets roughly that proportion of 
#		 1's in their lDel gene. Once their lDel gene is built, the rho value 
# 		 maximizing their fitness (see findRhoMaximizingFitness) is calculated. 
# 		 Then the individuals alpha and beta genes are initialized to all 0's. 
# 		 Finally, all of these are stored in a list and stored into the entry of 
#		 a larger list. Their index in the larger list is how they are referred to 
# 		 later on.
#
#	Arguments: population - the list where the genes and fitness of an individual are 
# 			        are stored
# 		   pOne - the proportion of 1's in the lDel gene for each individual in the 
#			  population
#		   s - an integer corresponding to the fitness cost associated with 
#		       proofreading
# 		   lDelGeneLength - the length of the lDel gene
#		   pNonDelToDel - the probability of a 1 changing to a 0 in lDel via mutation
#		   pDelToNonDel - the probability of a 0 changing to a 1 in lDel via mutation
#		   plDelLociMutation - the probability of a 1 or a 0 changing to a 0 or 1 
# 				       in lDel via mutation
# 		   alphaGeneLength - the number of floats in each alpha gene
#		   betaGeneLength - the number of floats in each beta gene
#
#	Returns: Nothing
#
##########################################################################################
def initializePopulation(populationSize, pOne, s, lDelGeneLength, 
                         pNonDelToDel, pDelToNonDel, plDelLociMutation, alphaGeneLength,
                         betaGeneLength):

	# Store the myPop as a list
	myPop = []                
	# Creates rho and lDel genes, in the form of two arrays, for each individual 
	# and calculates their fitness
	for i in range(populationSize):
		# Randomly lDel. Make approximately pOne of the lDels 1's and the rest 0's
		lDelGene = []
		for j in range(lDelGeneLength):
			if random.random() < pOneLoci:
				lDelGene.append(1)
			else:
				lDelGene.append(0)
	
		# Initialize alpha and beta genes to zeroes
		alphaGene = [0 for i in range(alphaGeneLength)]
		betaGene = [0 for i in range(betaGeneLength)]
		
		# Find the rho value maximizing fitness given the lDel of the individual
		individualRho = findRhoMaximizingFitness(s, lDelGene, pNonDelToDel, pDelToNonDel, 
		                                         plDelLociMutation, alphaGene, betaGene)
		
		# Calculate considering the lDel gene and rho gene 
		fitness = getFitness(s, lDelGene, pNonDelToDel, pDelToNonDel, plDelLociMutation, 
							 individualRho, alphaGene, betaGene, 0)
	
		# An individual consists of a rho gene, an lDel gene, and a fitness
		# Each is an entry in an myPopSize size array
		myPop.append([individualRho, lDelGene, fitness, alphaGene, betaGene])
	
	# Return a reference to the newly made population	
	return myPop
				
##########################################################################################
#	Name: pickDeadIndiv
#
#	Assumptions: It is assumed that list holding individuals, the ith entry in the list
#		     corresponds to the ith individual in the population. Furthermore, it is 
# 		     assumed the fitnessIndex entry in the ith index of the population array
# 		     is the fitness of the individual. Finally, it is assumed that the 
#		     fitness of an individual is less than one; if all indviduals in the 
#		     population have a fitness greater or equal to one, then an infinite 
#		     loop will occur.
#
#	Purpose: The simulation this function is used for is one that keeps a constant 
#		 population size. Therefore, to produce an offspring, an individual must
#		 "die" or be replaced. This function picks the dead indvidual. If selection
#		 is enabled, then a random index in the population representing an individual
#		 will be chosen. If selection is not on, then the first individual to have 
#		 a fitness less than a random number in the range [0, 1) will be picked as 
#		 the dead indvidual.
#
#	Arguments: selection - a boolean corresponding to whether or not selection will 
#			       dictate if the person who dies is random or someone whose
#			       fitness was not high enough
#		   population - a list containing the individuals of the population and all of
#				their genetic information
#		   fitnessIndex - the index in the an in individuals entry in the population 
#				  list holding the fitness of that individual. 
#				  I.e. population[indvidual][fitnessIndex] gives individual's
#				  fitness
#
#	Returns: An integer is returned. This integer is greater than zero and less than the
#		 number of individuals in the population. This integer corresponds to which
#		 individual has been chosen to die.
#
##########################################################################################
def pickDeadIndiv(selection, population, fitnessIndex):
	# Calculate how big the population is
	populationSize = len(population)
	
	# If selection is on, consider an individual's fitness
	if selection:
		# Pick a random individual until their fitness is lower than a random float 
		# between zero and one.
		while True:
			potientialDeadIndiv = int(random.random() * populationSize) 
			
			# If someone's fitness less than a random float, return their index in 
			# the population list
			if population[potientialDeadIndiv][fitnessIndex] < random.random():
				deadIndivIndex = potientialDeadIndiv
				break
		return deadIndivIndex
	# If selection is not on, then just pick a random person
	else: 
		return int(random.random() * populationSize) 
		
##########################################################################################
#	Name: replaceDeadWithOffspring
#
#	Assumptions: It is assumed that the likelihood of getting a rho gene from either
#		     in the case of recombination is .5; it is also assumed that sections
#		     of the lDel contributed from both parents can differ in size. In the case
#		     of no recombination, it is assumed that the person who died cannot be
#                    the parent of the offspring replacing them.
#
#	Purpose: This method exists so that a dead individual can be replaced with an 
#	         offspring. There are two cases: recombination is off or it is on. If it
#		 is off, then the offspring gets copies of each of the genes of another 
#		 individual in the population who is randomly chosen. If it is on, then 
#		 the offspring gets a combination of different genes from each parent and 
# 		 a fitness reflecting the difference in gene.
#			 
#
#	Arguments: deadIndex - the index in the population list of the individual who 
#			       died. This will be the index of the offspring.
#		   population - a list containing the individuals of the population and all of
#				their genetic information
#		   recombination - a boolean corresponding to whether or not recombination
#				   is allowed when developing offspring.
#
#	Returns: Nothing
#
##########################################################################################
def replaceDeadWithOffspring(deadIndex, recombination, population):
	if recombination:
		# Pick two parents
		mateOne = deadIndex
		while mateOne == deadIndex:
			mateOne = int(random.random() * populationSize)
		mateTwo = deadIndex
		while mateTwo == deadIndex:
			mateTwo = int(random.random() * populationSize)	
		
		# Randomly assign one of the rho values from the parents to the individual
		if random.random() < .5:
			population[deadIndex][0] = population[mateOne][0]
		else:
			population[deadIndex][0] = population[mateTwo][0]
		
		# Pick a recombination site for every fifty loci in the not including index first
		# or last index for the lDel gene
		numberOflDelRecombSites = int(lDelGeneLength / 50)
		lDelRecombSites = []
		lDelRecombSites.append(0)
		for i in range (numberOflDelRecombSites):
			lDelRecombSites.append(random.randint(0,lDelGeneLength - 1))
		lDelRecombSites.append(lDelGeneLength)
		lDelRecombSites.sort()
		
		# Splice lDel from recombination site to recombination site, alternating between
		# parent one and parent two
		population[deadIndex][1] = []
		for i in range (0, len(lDelRecombSites) - 1):
			if i % 2 == 1:
				population[deadIndex][1] += population[mateTwo][1]\
				                     [lDelRecombSites[i]:lDelRecombSites[i + 1]] 
			else:
				population[deadIndex][1] +=  population[mateOne][1]\
				                      [lDelRecombSites[i]:lDelRecombSites[i + 1]]
				                      
		# Pick recombination sites for alpha
		alphaLength = len(population[deadIndex][3])
		alphaRecombSites = 2 
		alphaDividers = []
		alphaDividers.append(0)
		for i in range (alphaRecombSites):
			alphaDividers.append(random.randint(0, alphaLength - 1))
		alphaDividers.append(alphaLength)
		alphaDividers.sort()
		
		# Splice together sections of alphas from both parents
		population[deadIndex][3] = []
		for i in range (0, len(alphaDividers) - 1):
			if i % 2 == 1:
				population[deadIndex][3] += population[mateTwo][3]\
				                            [alphaDividers[i]:alphaDividers[i + 1]] 
			else:
				population[deadIndex][3] += population[mateOne][3]\
											[alphaDividers[i]:alphaDividers[i + 1]]
											
		# Pick recombination sites for beta
		betaLength = len(population[deadIndex][4])
		betaRecombSites = 2
		betaDividers = []
		betaDividers.append(0)
		for i in range (betaRecombSites):
			betaDividers.append(random.randint(0, betaLength - 1))
		betaDividers.append(betaLength)
		betaDividers.sort()
		
		# Splice together beta genes
		population[deadIndex][4] = []
		for i in range (0, len(betaDividers) - 1):
			if i % 2 == 1:
				population[deadIndex][4] += population[mateTwo][4]\
				                            [betaDividers[i]:betaDividers[i + 1]] 
			else:
				population[deadIndex][4] +=  population[mateOne][4]\
								             [betaDividers[i]:betaDividers[i + 1]]
		
	else:
		# Pick an individual to be the parent who isn't the person who just died
		mateOne = deadIndex
		while mateOne == deadIndex:
			mateOne = int(random.random() * populationSize)
			
		# Assigne the give the offspring the rho of the parent
		population[deadIndex][0] = population[mateOne][0]
		
		# Copy each lDel loci from the parent to the lDel gene of the offspring.
		# (this is done manually to avoid multiple pointers to one gene)
		population[deadIndex][1] = []
		for i in range (len(population[mateOne][1])):
			population[deadIndex][1].append(population[mateOne][1][i])
		
		# Copy the alpha gene from the parent to the offspring
		population[deadIndex][3] = []
		for i in range(len(population[mateOne][3])):
			population[deadIndex][3].append(population[mateOne][3][i])
			
		# Copy the bea gene from the parent to the offspring
		population[deadIndex][4] = []
		for i in range(len(population[mateOne][4])):
			population[deadIndex][4].append(population[mateOne][4][i])
##########################################################################################
#	Name: mutateIndividual
#
#	Assumptions:
#
#	Purpose: This function exist for the sake of taking an individual and mutating their,
#                rho, lDel, alpha, and beta gene. A rho mutation takes an existing rho and 
#		 scales it. An lDel mutation changes a loci from a 1 to a 0 or a 0 to a 1
#		 (I don't know if it is possible to have no mutation). An alpha mutation 
#		 entails picking a random loci and incrementing or decrementing it. A beta
#		 mutation works the same as an alpha mutation except the the loci changed is
#	         in the beta gene. All of the mutations are done with rates given by the
#			  caller of the function. 
#
#	Arguments: population - the list storing all of the indviduals and their genetic 
#			        information
#		   plDelMutation - the probability of an
#		   pRhoMutation - the probability of a mutation in rho
# 		   pNonDelToDel - the probability of a 1 changing to a 0 in lDel via mutation
#		   pDelToNonDel - the probability of a 0 changing to a 1 in lDel via mutation
#		   mutantIndex - the index of the population  list where the genes of the 
#			         individual who will receive mutations. 
#				 population[mutantIndex] is the list of the mutants genes
#
#	Returns: None
#
##########################################################################################
def mutateIndividual(population, plDelMutation, pRhoMutation, pNonDelToDel, pDelToNonDel, 
                     mutantIndex):
	# Mutate rho with rate pRhoMutaion
	if random.random() < pRhoMutation:
		# I'm not too sure why rho is mutated this way
		population[mutantIndex][0] *= 10**random.gauss(0, .2)
	
	# Mutate lDel with rate pl	DelMutation
	if random.random() < plDelMutation:
		# The expression below after the less than sign is the probability of a zero 
		# going to a one, a loci going from non deleterious to deleterious
		if random.random() < population[mutantIndex][1].count(0) * pNonDelToDel / \
		                     (population[mutantIndex][1].count(0) * pNonDelToDel + \
		                     population[mutantIndex][1].count(1) * pDelToNonDel):# 0-->1
		    # The line below returns all of the indices in zeroes the offspring lDel 
		    # gene
			lDelZeroIndices = [i for i, x in enumerate(population[mutantIndex][1])\
			                   if x == 0]
			# This condition is here to check for the case that there's a mutation from
			# a zero to a one but no zeroes exists in the lDel gene
			if len(lDelZeroIndices) > 0:
				# The line below changes a random 0 to a one
				population[mutantIndex][1][lDelZeroIndices[int(random.random() * \
				                                  len(lDelZeroIndices))]] = 1
		else: # That mutation 1->0
			#print("PURGED\n")
			# The line below returns all of the indices in ones the offspring lDel 
		    # gene
			lDelOneIndices = [i for i, x in enumerate(population[mutantIndex][1]) \
			                  if x == 1]
			# This condition is here to check for the case that there's a mutation from
			# a one to a zero but no ones exists in the lDel gene
			if len(lDelOneIndices) > 0:
				# The line below changes a random 1 to a 0
				population[mutantIndex][1][lDelOneIndices[int(random.random() * \
				                  len(lDelOneIndices))]] = 0
				                  
		'''
		if random.random() < imNotSureYet:
			newAlpha[random.random() * len(newAlpha)] += \
			random.gauss((-1 * sum(newAlpha) / 50), 1)
		if random.random() < someOtherNumber	
			newBeta[random.random() * len(newBeta)] += \
			random.gauss((-1 * sum(newBeta) / 50), 1)
		'''
##########################################################################################
#	Name: outputlDelCount
#
#	Assumptions: It is assumed that the individuals in a population are stored in a 
#		     list. Furthermore, it is assumed that each entry in the population
#		     list is another list containing genes of a given individual. 
#	             Finally, it is assumed that an lDel gene is a list and each 1 in the
#		     lDel gene is a loci that is deleterious.
#
#	Purpose: The purpose of this method is to output the lDel count to a file. This
#		 file will then be used to make a density plot in R. To accomplish this,
#		 each individual's lDel gene will be counted for lDels and then written
#		 to the desired output file on a new line. To be able to distinguish 
#		 between lDels of on group and another in the output file, a blank line
#		 will be used as a divider in the file, hence the new line write at the
#		 end of the function.
#
#	Arguments: - population is a list containing the lDels of each idividual
#		   - lDelIndex is the index of the sub list where the lDel of an indivudal
#		     can be accessed. I.e. population[individual][lDelIndex] gives access
#		     to individual's lDel gene
#		   - outfile is the file object representing the text file that the lDel
#		     measurements will be written to.
#
#	Returns: Nothing
#
##########################################################################################
def outputlDelCount(population, lDelIndex, outfile):
	for i in range (len(population)):
		outfile.write(str(population[i][lDelIndex].count(1)) + " ")
	# The line below ensures that each group of ldel outputs is 
	# seperated by a new line
	outfile.write("\n")
#-----------------------------END FUNCTION DEFINITIONS------------------------------------
# Output the time the program started
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Program started...\n")
sys.stdout.flush()

# All of the parameters that are set via the command line
recombination = int(sys.argv[1])
os.chdir(sys.argv[2])
populationSize = int(sys.argv[3])
lDelGeneLength = int(sys.argv[4])
generations = int(sys.argv[5])
pOneLoci = float(sys.argv[6])
pZeroLoci = 1 - pOneLoci

# DEBUGGING
print("Recomb", recombination, "Directory", sys.argv[2], "PopSize", populationSize, "# of lDels", lDelGeneLength, "Gens", generations, "pOne", pOneLoci)
# DEBUGGING

# All of the other paramaters. Most of these are taken from the 2011 paper
selection = True
pNonDelToDel = .4
pDelToNonDel = .1 
# per base mutation rate
plDelLociMutation = 10**-8
delLociFitnessCost = 7
# multiply below times 1000
pRhoMutation = 10**-5
nucleotidesPerlDel = 30
plDelMutation = plDelLociMutation * nucleotidesPerlDel * lDelGeneLength
alphaGeneLength = 10
betaGeneLength = alphaGeneLength
pAlphaMutation = plDelLociMutation * nucleotidesPerlDel * alphaGeneLength
pBetaMutation = plDelLociMutation * nucleotidesPerlDel * betaGeneLength
envOpt = 0

# How man times statstics are to be recorded to text files
NUMBER_OF_RESULT_WRITES = 100
NUMBER_OF_LDEL_COUNT_WRITES = 5

# Call for the population to be initialized
population = initializePopulation(populationSize, pOneLoci, delLociFitnessCost, 
								  lDelGeneLength, pNonDelToDel, pDelToNonDel, 
								  plDelLociMutation, alphaGeneLength, betaGeneLength)
'''								  
# DEBUGGING
print("Initial Population")
for i in range (len(population)):
	print(i, ": ", population[i][0], len(population[i][1]), population[i][2], "\n")
# DEBUGGING
'''
# Output the time the population is done being initialized 
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Population initialized...\n")
sys.stdout.flush() 
# Open and write the columns for the results text file
results = open("results.txt", "w")
columns = ["MeanFitness", "MeanFitnessVariance", "MeanlDels", "MeanlDelsVariance",
           "MeanRho", "MeanRhoVariance", "MeanAlpha", "MeanAlphaVariance", 
           "MeanBeta", "MeanBetaVariance"]
results.write("\t".join(columns) + "\n")
lDelOut = open("lDelOutput.txt", "w")
	
# Pick an individual to die, produce an offspring, mutates the offspring and report
# statstics on the population. One generation is N deaths and replacements
for replacementNumber in range (populationSize * generations):
	deadIndex = pickDeadIndiv(selection, population, 2)	
	'''
	# DEBUGGING
	if replacementNumber % (populationSize * generations / 100) == 0:
		print("=============================================================================")
		print("Dead Individual: ", deadIndex, population[deadIndex][0], len(population[deadIndex][1]), population[deadIndex][2], "\n")
	# DEBUGGING
	'''	
	replaceDeadWithOffspring(deadIndex, recombination, population)
	'''
	# DEBUGGING
	if replacementNumber % (populationSize * generations / 100) == 0:
		print("New Individual: ", deadIndex, population[deadIndex][0], len(population[deadIndex][1]), population[deadIndex][2], "\n")
	# DEBUGGING
	'''
	mutateIndividual(population, plDelMutation, pRhoMutation, pNonDelToDel, pDelToNonDel,
	                 deadIndex)
	'''
	# DEUBGGING
	if replacementNumber % (populationSize * generations / 100) == 0:
		print("New Mutated Individual: ", deadIndex, population[deadIndex][0], len(population[deadIndex][1]), population[deadIndex][2], "\n")
	# DEBUGGING
	'''
	population[deadIndex][2] = getFitness(delLociFitnessCost, population[deadIndex][1], 
	                                      pNonDelToDel, pDelToNonDel, plDelMutation, 
	                                      population[deadIndex][0], 
	                                      population[deadIndex][3], 
	                                      population[deadIndex][4], envOpt)
	'''                                      
	# DEBUGGING
	if replacementNumber % (populationSize * generations / 100) == 0:
		print("Paramaters of fitness function")
		print ("delLociFitnessCost", delLociFitnessCost)
		print("lDel", population[deadIndex][1])
		print("pDel", pNonDelToDel)
		print("pMinsDel", pDelToNonDel)
		print("plDelMutation", plDelMutation)
		print("rho", population[deadIndex][0])
		print("alpha", population[deadIndex][3])
		print("beta", population[deadIndex][4])
		print("Environmental Optimum", envOpt)
		print("Fitness of offspring", population[deadIndex][2])
		print("=============================================================================")
	# DEBUGGING
	'''
	
	# Records stats 100 times every run and report the progress of the script
	if replacementNumber % (populationSize * generations / NUMBER_OF_RESULT_WRITES) == 0:
		recordStatistics(results, population)

		sys.stdout.write("[" + time.asctime() + "]: ")
		sys.stdout.write(str(int(replacementNumber / \
		                 (populationSize * generations) * NUMBER_OF_RESULT_WRITES)) + \
		                 "% complete\n")
		sys.stdout.flush()
		
	if replacementNumber % (populationSize * generations / \
	   NUMBER_OF_LDEL_COUNT_WRITES) == 0:
		outputlDelCount(population, 1, lDelOut)

# Write the parameters to the end of the text file and close the file
parameters = "[populationSize = " +  str(populationSize) +  ", Generations = " + \
             str(generations) +  ", pOneLoci = " +  str(pOneLoci) + \
             ", Recombination = "  + str(recombination) + ", Ldel length = " + \
             str(lDelGeneLength) +"]" 
results.write(parameters)
results.write("\n")
results.close()
lDelOut.close();

# Pickle the population into a binary file
with open("frozenPopulation.bin", "wb") as storage:
	pickle.dump(population, storage)

# Output the time that the program finishes
sys.stdout.write("[" + time.asctime() + "]: ")
sys.stdout.write("Program done\n")
sys.stdout.flush()#aasdfgasfgads
