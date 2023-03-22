#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rng.h"		
#include "rng.c"

// NOTES FOR GRADING:
// I chose to implement intermediate recombination since I felt that this kind of function optimization
// problem could benefit from some extra exploitation. I tested with and without a non-zero crossover
// rate, however, and didn't notice any significant improvement in convergence rate. Implementation
// of the actual crossover is in "ES_Genome_Recombo(...)" and it gets used in
// "ES_Population_Make_New_Generation(...)" when CROSSOVER_RATE > 0.

// I also chose to use (mu+lambda) survivor selection since parent fitness may often be better than child
// fitness for this kind of stochastic hill climber, especially in early generations before appropriate
// sigma values have been selected. I used rank order survivor selection with the exponential probability
// formula given by the textbook to apply more selection pressure than a simple linear distribution. I did
// modify the given formula by effectively scaling down the input (index) exponent by the size of the
// population being selected from, since the curve of the given probability function leveled out way too 
// quickly to be useful for a selection pool of over a hundred genomes. This can all be seen in 
// "ES_Population_Make_New_Generation(...)".

// I used correlated mutations (individual sigmas for each allele) without rotation, frankly because it was
// the easier option. I also feel like this simple of an optimization problem doesn't warrant the extra
// performance/calculation cost of rotation parameters, which I think is evident based on how quickly this
// version converges. Mutation implementation is in "ES_Genome_Mutate(...)" and can also be seen implemented
// separately in "ES_Genome_Init(...)".

// ------- Global EC Parameters -------

// const int LAMBDA_SIZE       = 100;      // Size of child "population" (# of offspring per generation)
// const int MU_SIZE           = 15;       // Size of ongoing/maintained parent population
// const int GENOME_LEN        = 2;        // # of floating pt genes *NOT including sigmas*
// const double MUTATION_RATE  = 0.5;      // Mutation rate parameter for sigmas ("Tau")
// const double CROSSOVER_RATE = 0.0;
// const int GENERATION_COUNT  = 1000;     // Max # of generational cycles

// // NOTE: for simplicity I'm using the same range constants for each (both, in this case) genes
// const double RANGE_MIN      = -10.0;    // Min value of each gene
// const double RANGE_MAX      = 10.0;     // Min value of each gene

// // This parameter keeps any sigma from getting so close to 0 that there's effectively no mutation
// const double SIGMA_MIN      = 0.00001;


// ------- Data Structures and Type Definitions -------

// TODO: change to proper bit string representation
typedef struct genome_s {
    int genome_len; 	    // genome length (size of array of doubles)
    double *gene;           // pointer to an array of doubles          
	double fitness;         // fitness score of genome
} genome_t;


typedef struct population_s {
    int member_count;       // number of members (genomes) in pop.
	genome_t *member;       // An array of pointers to genome structures
	} population_t;
		
// Function prototype to keep GCC happy with the order of definitions/calls
void ES_Population_Compute_Fitness(population_t *);


// ------- Population Memory Allocation Routines -------

// allocate memory for a genome of <length> doubles
void ES_Genome_Malloc(genome_t *genome, int length) {
    (*genome).gene = (double *)malloc(sizeof(double)*length);
    (*genome).genome_len = length;
    (*genome).fitness = 0.0;
}

//deallocate memory used by a genome
void ES_Genome_Free(genome_t *genome) {
    free((*genome).gene);
    (*genome).genome_len = 0;
}

// allocate population of <member_count> floating pt genomes of length <genome_length>
void ES_Population_Malloc(population_t *population, int member_count, int genome_length) {
    int member_c;

    (*population).member_count = member_count;
    (*population).member = (genome_t *)malloc(member_count * sizeof(genome_t));
    for (member_c=0; member_c < member_count; member_c++)
	    ES_Genome_Malloc(&((*population).member[member_c]), genome_length);	
}


void ES_Population_Free(population_t *population) {
    // modeled after the provided functions to deallocate all malloc'd genomes,
    // then deallocate the array of genome pointers itself
    for (int member_c=0; member_c < (*population).member_count; member_c++)
	    ES_Genome_Free(&((*population).member[member_c]));
    
    free((*population).member);
    (*population).member_count = 0;
}


// ------- Population Utility Routines -------

// print out genome in brackets followed by fitness score
void ES_Genome_Print(genome_t *genome){
    int gene_count;
    
    printf("[");
    for (gene_count = 0; gene_count < (*genome).genome_len - 1; gene_count++) {
        printf("%f,", (*genome).gene[gene_count]);
    }
    printf("%f] %f\n", (*genome).gene[gene_count], (*genome).fitness);
    
	
}

// Print ALL the genomes in a population, with fitness scores
void ES_Population_Print(population_t *population) {
    int p_count;
    
    for (p_count=0; p_count < (*population).member_count; p_count++) {
        ES_Genome_Print(&((*population).member[p_count]));
    }
    printf("\n");
}

// Copy info in <*source_genome> to <*destination_genome>
void ES_Genome_Copy(genome_t *source_genome, genome_t *destination_genome) {
    // NOTE: assumes that both have allocated memory
    // TODO: add error detection for memory allocation?
    int gene_pos;
    if ((*source_genome).genome_len != (*destination_genome).genome_len) {
        fprintf(stderr, "ERROR in ES_Genome_Copy: Source and Destination of Unequal Length\n");
        exit(0);
    }
 
    (*destination_genome).fitness = (*source_genome).fitness;
    for (gene_pos=0; gene_pos<(*destination_genome).genome_len; gene_pos++) {
        (*destination_genome).gene[gene_pos] = (*source_genome).gene[gene_pos];
    }
}

// Copy source population to dest. population (source unmodified)
void ES_Population_Copy(population_t *source_population, population_t *destination_population) {
    // TODO: add error detection for genome size within populations
    int p_count;
    if ((*source_population).member_count != (*destination_population).member_count) {
        fprintf(stderr, "ERROR in ES_Population_Copy: Source and Destination of Unequal Sizes\n");
        exit(0);
    }
    for (p_count=0; p_count < (*source_population).member_count; p_count++) {
        ES_Genome_Copy(&((*source_population).member[p_count]),&((*destination_population).member[p_count]));
    }
}

// Calculate metrics (Max fitness, avg. fitness, diversity score, champ x/y)
void ES_Population_Metrics(population_t *population, double *metrics) {
    int m_count;
    int m_inner;
    int g_bit;
    double sum_of_fitness = 0.0;
    int num_unique = 0;
    int champ_id;
    double x1 = 0.0;
    double y1 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    double distance;

    // move previous generation's metrics to end of array
    // (prev. metrics kept for possible rate convergence calculations - not currently used)
    metrics[5] = metrics[0];
    metrics[6] = metrics[1];
    metrics[7] = metrics[2];
    metrics[8] = metrics[3];
    metrics[9] = metrics[4];

    // find max fitness score and sum of all scores for avg calc
    double max_fitness = ((*population).member[0]).fitness;
    champ_id = 0;
	for (m_count = 0; m_count < (*population).member_count; m_count++) {
        sum_of_fitness += ((*population).member[m_count]).fitness;

        if (((*population).member[m_count]).fitness > max_fitness) {
            max_fitness = ((*population).member[m_count]).fitness;
            champ_id = m_count;
        }
    }	    

    // store max score in metrics
    metrics[0] = max_fitness;

    // store avg fitness of generation in metrics
    metrics[1] = (sum_of_fitness / (float)(*population).member_count);

    // store x, y of champion
    metrics[3] = (*population).member[champ_id].gene[0];
    metrics[4] = (*population).member[champ_id].gene[1];

    // Calculate furthest distance between members as diversity metric
    // DISABLED FOR PERFORMANCE WHILE PARAMETER TUNING
    metrics[2] = 0.0;
    // for (m_count = 0; m_count < (*population).member_count; m_count++) {
    //     x1 = (*population).member[m_count].gene[0];
    //     y1 = (*population).member[m_count].gene[1];
        
    //     for (m_inner = m_count+1; m_inner < (*population).member_count; m_inner++) {
    //         x2 = (*population).member[m_inner].gene[0];
    //         y2 = (*population).member[m_inner].gene[1];
    //         distance = sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0));

    //         // Store largest distance
    //         if (distance > metrics[2]) metrics[2] = distance;
    //     }
    // }
}

// qsort comparison function for population sorting by fitness
int ES_Qsort_Comp(const void *a, const void *b) {
    if(((genome_t*)a)->fitness > ((genome_t*)b)->fitness) return 1;
    if(((genome_t*)a)->fitness < ((genome_t*)b)->fitness) return -1;
    return 0;
}

// ------- Population Initialization Routines -------

// set *already allocated* genome to random set of doubles in search range
void ES_Genome_Init(genome_t *genome, double range_min, double range_max, double mutation_rate, RNG *rng) {
    int gene_count;
    for (gene_count=0; gene_count < (*genome).genome_len; gene_count++) {
        // set each allele to uniform random val between RANGE_MIN and RANGE_MAX (random resetting)
        if(gene_count < (*genome).genome_len / 2) {
            (*genome).gene[gene_count] = rng_uniform01(rng) * (range_max - range_min) + range_min;
        }
        // set each sigma to random val with lognormal distribution
        else {
            (*genome).gene[gene_count] = exp(mutation_rate * rng_gaussian(rng, 0, 1));    // textbook formula 4.2
        }
    }
    (*genome).fitness = 0.0;
}

// randomize all genomes in population
void ES_Population_Init(population_t *population, double range_min, double range_max, double mutation_rate, RNG *rng) {
    int m_count;
    for (m_count=0; m_count < (*population).member_count; m_count++) {
        ES_Genome_Init(&((*population).member[m_count]), range_min, range_max, mutation_rate, rng);
    }
}


// ------- Genome Decoding Functions -------

// NOTE: Since we're using floating pt representation for a simple function,
// our phenotype and genotype are the same thing and no decoding is needed


// ------- Genome Variation Functions -------

// mutate each sigma, then mutate each gene based on respective sigma value
void ES_Genome_Mutate(genome_t *genome, double mutation_rate, double range_min, double range_max, double sigma_min, RNG *rng) {
    int gene_count;
    double sigma_prime;
    double sigma_i;
    // Loop through genome backwards to modify sigmas first
    for (gene_count = (*genome).genome_len - 1; gene_count >= 0; gene_count--) {
        if (gene_count >= (*genome).genome_len / 2) {
            // NOTE: I've chosen here to use the "simpler" equation for mutating
            // each sigma based on a single Tau (i.e. textbook equation 4.2) 
            // rather than the more "granular" equation 4.4.
            sigma_prime = (*genome).gene[gene_count] * exp(mutation_rate * rng_gaussian(rng, 0, 1));
            if (sigma_prime < sigma_min)
                sigma_prime = sigma_min;
            (*genome).gene[gene_count] = sigma_prime;
        }
        // Now mutate each allele based on its individual (already mutated) sigma
        else {
            sigma_i = (*genome).gene[gene_count + (*genome).genome_len / 2];
            (*genome).gene[gene_count] += sigma_i * rng_gaussian(rng, 0, 1);

            // Bound genes to defined range
            if ((*genome).gene[gene_count] > range_max) (*genome).gene[gene_count] = range_max;
            else if ((*genome).gene[gene_count] < range_min) (*genome).gene[gene_count] = range_min;
        }
    } 
}

// Recombination of two genomes - intermediate
void ES_Genome_Recombo(genome_t *parent_genome_1, genome_t *parent_genome_2, RNG *rng) {
    // NOTE: I coded this up just to think through it, but my crossover_rate is 0
    // so this isn't currently getting used
    int gene_count;
    double alpha;
    genome_t child_genome_1, child_genome_2;

    // malloc children genomes with same length as parents
    ES_Genome_Malloc(&child_genome_1, (*parent_genome_1).genome_len);
    ES_Genome_Malloc(&child_genome_2, (*parent_genome_1).genome_len);

    // apply intermediate recombination to kids with uniform random alpha for alleles and sigmas
    for (gene_count = 0; gene_count < (*parent_genome_1).genome_len; gene_count++) {
        alpha = rng_uniform01(rng);
        (child_genome_1).gene[gene_count] = alpha * (*parent_genome_1).gene[gene_count] + (1 - alpha) * (*parent_genome_2).gene[gene_count];
        (child_genome_2).gene[gene_count] = (1 - alpha) * (*parent_genome_1).gene[gene_count] + alpha * (*parent_genome_2).gene[gene_count];
    }
	
    // replace parents with children
    ES_Genome_Copy(&child_genome_1, parent_genome_1);
    ES_Genome_Copy(&child_genome_2, parent_genome_2); 

    // deallocate temp child memory
    ES_Genome_Free(&child_genome_1);
    ES_Genome_Free(&child_genome_2); 
}

// apply mutation to every genome in population
void ES_Population_Mutate(population_t *population, double mutation_rate, double range_min, double range_max, double sigma_min, RNG *rng) {
    // TODO: add error checking?
    int m_count;
	
    // loop through every member of population and apply ES_Genome_Mutate
    for (m_count = 0; m_count < (*population).member_count; m_count++) {
        ES_Genome_Mutate(&((*population).member[m_count]), mutation_rate, range_min, range_max, sigma_min, rng);
    }
}

void ES_Population_Make_New_Generation(population_t     *old_population, 
                                        population_t    *new_population,
                                        int             genome_len,
                                        double          mutation_rate,
                                        double          crossover_rate,
                                        double          mu_size,
                                        double          lambda_size,
                                        double          range_min,
                                        double          range_max,
                                        double          sigma_min,
                                        double          *metrics,
										RNG             *rng) {

    double sum_of_fitness_scores;
    double max_raw_fitness, min_raw_fitness;

    int m_count;
    int parent_one_index;
    int parent_two_index;
    double *rank;
    double norm_factor_c;
    double die_roll;

    population_t selection_pool;

    // Malloc temp population to hold Mu + Lambda
    ES_Population_Malloc(&selection_pool, (mu_size + lambda_size), (genome_len * 2));

    // Malloc array for rank order PDF/CDF construction
    rank = (double *)malloc(sizeof(double)*(mu_size + lambda_size));

    // Copy parent population into selection pool
    for (m_count = 0; m_count < mu_size; m_count++) {
        ES_Genome_Copy(&((*old_population).member[m_count]), &(selection_pool.member[m_count]));
    }

    // Generate Lambda children from Mu parents using uniform random selection, fill rest of pool

    // [NOTE: We generate kids 2 at a time to enable possible crossover, but the way we do this means
    // that for every loop iteration but the first/last we overwrite a child we already generated.
    // This is pretty memory inefficient but allows the use of odd population sizes without a bunch
    // of conditionals.]
    for (; m_count < (mu_size + lambda_size)-1; m_count++) {
        // Uniform random selection of two parents
        parent_one_index = (int)floor(rng_uniform(rng, 0.0, (double)(*old_population).member_count));
        parent_two_index = (int)floor(rng_uniform(rng, 0.0, (double)(*old_population).member_count));

        // Clone parents to generate kids before (possible) mutation
        ES_Genome_Copy(&((*old_population).member[parent_one_index]), &(selection_pool.member[m_count]));
        ES_Genome_Copy(&((*old_population).member[parent_two_index]), &(selection_pool.member[m_count+1]));

        // Mutate the children
        ES_Genome_Mutate(&(selection_pool.member[m_count]), mutation_rate, range_min, range_max, sigma_min, rng);
        ES_Genome_Mutate(&(selection_pool.member[m_count+1]), mutation_rate, range_min, range_max, sigma_min, rng);

        // Recombine the two generated children
        if(rng_uniform01(rng) < crossover_rate) {
            ES_Genome_Recombo(&(selection_pool.member[m_count]), &(selection_pool.member[m_count+1]), rng);
        }
    }

    // ~~~ RANK ORDER, MU+LAMBDA SURVIVOR SELECTION ~~~
    
    // Update fitness of selection_pool members
    ES_Population_Compute_Fitness(&selection_pool);

    // Sort selection_pool.member[] via qsort and fitness comparison function
    qsort(selection_pool.member, selection_pool.member_count, sizeof(genome_t), ES_Qsort_Comp);

    // Calculate normalization factor "c" 
    norm_factor_c = 0.0;
    for (m_count = 0; m_count < selection_pool.member_count; m_count++) {
        norm_factor_c += (1.0 - exp(-1.0*m_count/(double)selection_pool.member_count));
    }

    // Use c to construct discrete PDF in rank[]
    for (m_count = 0; m_count < selection_pool.member_count; m_count++) {
        rank[m_count] = (1.0 - exp(-1.0*m_count/(double)selection_pool.member_count)) / norm_factor_c;
    }

    // Construct CDF by adding up probabilities
    for (m_count = 1; m_count < selection_pool.member_count; m_count++) {
        rank[m_count] += rank[m_count-1];
    }

    // Select survivors placed in new_population via CDF
    for (m_count = 0; m_count < (*new_population).member_count; m_count++) {
        die_roll = rng_uniform01(rng);
        parent_one_index = 0;
        while (rank[parent_one_index] < die_roll) parent_one_index++;
        ES_Genome_Copy(&(selection_pool.member[parent_one_index]), &((*new_population).member[m_count]));
    }

    // Free memory
    ES_Population_Free(&selection_pool);
    free(rank);
}


// ------- Genome Fitness Functions -------

double ES_Rosenbrock(genome_t *genome) {
    const double a = 1.0;
    const double b = 100.0;
    
    double x = (*genome).gene[0];
    double y = (*genome).gene[1];

    double f_xy = pow((a - x), 2.0) + b * pow((y - x*x), 2.0);

    return 0.0 - f_xy;
}

double ES_Himmelblau(genome_t *genome) {
    double x = (*genome).gene[0];
    double y = (*genome).gene[1];

    double f_xy = pow(x*x - 11 + y, 2.0) + pow((x + y*y - 7), 2.0);
    return (0.0 - f_xy);
}

// wrapper for desired fitness function
double ES_Fitness_Function(genome_t *genome) {
    return ES_Rosenbrock(genome);
}

// apply fitness function wrapped by ES_Fitness_Function to whole population
void ES_Population_Compute_Fitness(population_t *population) {
    int m_count;
	
    for (m_count = 0; m_count < (*population).member_count; m_count++)
        ((*population).member[m_count]).fitness = ES_Fitness_Function(&((*population).member[m_count]));

}

// -----------------------------------------------

int main() {

    // ------- Global EC Parameters -------

    int lambda_size         = 100;      // Size of child "population" (# of offspring per generation)
    int mu_size             = 15;       // Size of ongoing/maintained parent population
    int genome_len          = 2;        // # of floating pt genes *NOT including sigmas*
    double mutation_rate    = 0.5;      // Mutation rate parameter for sigmas ("Tau")
    double crossover_rate   = 0.0;

    double range_min      = -5.11;    // Min value of each gene
    double range_max      = 5.12;     // Min value of each gene

    // This parameter keeps any sigma from getting so close to 0 that there's effectively no mutation
    double sigma_min      = 0.001;

    // ------- Test Parameters -------

    int generation_count    = 1000;     // Max # of generational cycles
    double optimum_fitness  = -0.001;
    int run_count           = 100;

    double mutation_min     = 0.2;
    double mutation_max     = 0.6;
    double sigma_ctrl_min   = 0.001;
    double sigma_ctrl_max   = 0.02;
    double crossover_min    = 0.0;
    double crossover_max    = 0.4;
    int step_count          = 5;        // Number of values tested for each parameter INCLUDING end points

    // -------
    
    RNG *random_number_generator;
	population_t POPULATION,POPULATION2;
	int generation;
    double generation_avg = 0.0;    // Stores running sum of generations to reach optimum (or not) for avg calculation
    double champ_avg = 0.0;         // ^^^ Same for best fitness
    double champ_stddev = 0.0;
    FILE *csv;

    // allocate memory for array to hold running test metrics:
    // max fitness, avg fitness, diversity score, champ x, champ y (for generation n and n-1)
    double metrics[10] = {0.0};
	
	// initialize jcg's rng
	random_number_generator = rng_create();
	
	// allocate memory for populations with defined sizes
	ES_Population_Malloc(&POPULATION, mu_size, genome_len*2);  // primary population (GENOME_LEN doubled to account for sigmas)
	ES_Population_Malloc(&POPULATION2, mu_size, genome_len*2); // "scratch" population 
	
    for (mutation_rate = mutation_min; mutation_rate < (mutation_max*1.0001); mutation_rate += ((mutation_max - mutation_min)/(double)(step_count-1))) {
        printf("[mutation_rate = %f]\n", mutation_rate);

        for (crossover_rate = crossover_min; crossover_rate < (crossover_max*1.0001); crossover_rate += ((crossover_max - crossover_min)/(double)(step_count-1))) {
            printf("\t\t[crossover_rate = %f]\n", crossover_rate);

            for (sigma_min = sigma_ctrl_min; sigma_min < (sigma_ctrl_max*1.0001); sigma_min += ((sigma_ctrl_max - sigma_ctrl_min)/(double)(step_count-1))) {
                printf("\t\t\t[sigma_min = %f]\n", sigma_min);

                generation_avg = 0.0;
                champ_stddev = 0.0;
                champ_avg = 0.0;

                for (int i = 0; i < run_count; i++) {
                    // initialize/randomize population, calculate initial fitness and metrics
                    ES_Population_Init(&POPULATION, range_min, range_max, mutation_rate, random_number_generator);
                    
                    // ES_Population_Compute_Fitness(&POPULATION);
                    // ES_Population_Metrics(&POPULATION, metrics);
                    // // printf("Initial Population\n");
                    // ES_Population_Print(&POPULATION);

                    // generation evolutionary loop
                    for (generation = 0; generation < generation_count; generation++) {
                        // Evaluate Everyone in the population
                        ES_Population_Compute_Fitness(&POPULATION);

                        // Calculate and print metrics - fitness score of champ, avg fitness score, % identical genomes
                        ES_Population_Metrics(&POPULATION, metrics);
                        // printf("Rosenbrock ES %d %d %f %f %d %d %f %f %f\n", mu_size, lambda_size, mutation_rate, crossover_rate,
                        //                                                     generation, generation*(mu_size+lambda_size), 
                        //                                                     metrics[0], metrics[1], metrics[2]);
                        
                        // ES_Population_Print(&POPULATION);

                        // // Check for convergence - near-optimal solution
                        // if (metrics[0] > -0.000001) break;
                        
                        // Make the new generation
                        ES_Population_Make_New_Generation(  &POPULATION,
                                                            &POPULATION2, 
                                                            genome_len, 
                                                            mutation_rate, 
                                                            crossover_rate, 
                                                            mu_size,
                                                            lambda_size,
                                                            range_min,
                                                            range_max,
                                                            sigma_min,
                                                            metrics, 
                                                            random_number_generator);
                        // Copy the new population back into the primary population
                        ES_Population_Copy(&POPULATION2, &POPULATION);
                    }

                    if (generation == generation_count) {
                        // Calculate and print metrics for final generation if no convergence
                        ES_Population_Compute_Fitness(&POPULATION);
                        ES_Population_Metrics(&POPULATION, metrics);
                        // printf("Rosenbrock ES %d %d %f %f %d %d %f %f %f\n", mu_size, lambda_size, mutation_rate, crossover_rate,
                        //                                                     generation_count, generation_count*(mu_size+lambda_size), 
                        //                                                     metrics[0], metrics[1], metrics[2]);
                    }
                    
                    // printf("\n");	

                    generation_avg += generation;
                    champ_stddev += pow(metrics[0], 2.0);
                    champ_avg += metrics[0];

                    printf("\t\t\t[%f, %f]: %f; gen avg = %f\n", metrics[3], metrics[4], metrics[0], metrics[1]);
                }

                generation_avg /= run_count;
                champ_avg /= run_count;
                champ_stddev = sqrt((champ_stddev / (double)run_count) - pow(champ_avg, 2.0));

                csv = fopen("es_results.csv", "a");

                if (csv == NULL) {
                    printf("Error opening sga_results.csv");
                    exit(1);
                }
                fflush(csv);
                fprintf(csv, "%f, %f, %f, %f, %f, %f\n", generation_avg, mutation_rate, crossover_rate, sigma_min, champ_avg, champ_stddev);
                fclose(csv);
            
            }
        }
    }

	// printf("Final Population\n");
    // ES_Population_Print(&POPULATION); 
	
    // Print reason for stop
    // if (generation_count < GENERATION_COUNT)
    //     printf("\n[Terminated due to convergence]");
    // else
    //     printf("\n[Termination due to max generation count reached (%d)]", GENERATION_COUNT);
    
    return 0;
}
