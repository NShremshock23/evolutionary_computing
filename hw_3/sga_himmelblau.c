#include <stdio.h>
#include <stdlib.h>
#include "rng.h"		
#include "rng.c"

// NOTE FOR GRADING:
// All implementations are commented (in caps) with their applicable
// name, so you can easily find each by CTRL+F'ing
// (ex. the function containing roulette wheel is commmented with
// "ROULETTE WHEEL IMPLEMENTATION")

// Descriptive comments and section divisions similar to those in the
// demo code are also included.


// ------- Global EC Parameters -------

const int POP_SIZE          = 100;
const int GENOME_LEN        = 20;
const double MUTATION_RATE  = 0.002;
const double CROSSOVER_RATE = 0.9;
const int GENERATION_COUNT  = 1000;

const double RANGE_MIN      = -10.0;
const double RANGE_MAX      = 10.0;


// ------- Data Structures and Type Definitions -------

typedef struct genome_s {
    int genome_len; 	    // genome length in bits
	char *bit;              // pointer to an array of chars          
	double fitness;         // fitness score of genome
} genome_t;


typedef struct population_s {
    int member_count;       // number of members (genomes) in pop.
	genome_t *member;       // An array of pointers to genome structures
	} population_t;
		
// function prototype to keep GCC happy with the ordering of definitions and calls
double SGA_Genome_Decode(genome_t *, int, int, double, double);

// ------- Population Memory Allocation Routines -------

// allocate memory for a genome of <length> bits
void SGA_Genome_Malloc(genome_t *genome, int length) {
    (*genome).bit = (char *)malloc(sizeof(char)*length);
    (*genome).genome_len = length;
    (*genome).fitness = 0.0;
}

//deallocate memory used by a genome
void SGA_Genome_Free(genome_t *genome) {
    free((*genome).bit);
    (*genome).genome_len = 0;
}

// allocate population of <member_count> bit string genomes of length <genome_length>
void SGA_Population_Malloc(population_t *population, int member_count, int genome_length) {
    int member_c;

    (*population).member_count = member_count;
    (*population).member = (genome_t *)malloc(member_count * sizeof(genome_t));
    for (member_c=0; member_c < member_count; member_c++)
	    SGA_Genome_Malloc(&((*population).member[member_c]), genome_length);	
}


void SGA_Population_Free(population_t *population) {
    // modeled after the provided functions to deallocate all malloc'd genomes,
    // then deallocate the array of genome pointers itself
    for (int member_c=0; member_c < (*population).member_count; member_c++)
	    SGA_Genome_Free(&((*population).member[member_c]));
    
    free((*population).member);
    (*population).member_count = 0;
}


// ------- Population Utility Routines -------

// print out genome as a bitstring with current fitness score in parenthesis
void SGA_Genome_Print(genome_t *genome){
    int bit_count;
    double x = SGA_Genome_Decode(genome, 0, (genome->genome_len/2)-1, RANGE_MIN, RANGE_MAX);
    double y = SGA_Genome_Decode(genome, (genome->genome_len/2), genome->genome_len-1, RANGE_MIN, RANGE_MAX);

    for (bit_count=0; bit_count < (*genome).genome_len; bit_count++) {
        printf("%d",(*genome).bit[bit_count]);
    }
    printf(" [%f, %f] (%lf)\n", x, y, (*genome).fitness);
	
}

// Print ALL the genomes in a population, with fitness scores
void SGA_Population_Print(population_t *population) {
    int p_count;
    
    for (p_count=0; p_count < (*population).member_count; p_count++) {
        SGA_Genome_Print(&((*population).member[p_count]));
    }
}

// Copy info in <*source_genome> to <*destination_genome>
void SGA_Genome_Copy(genome_t *source_genome, genome_t *destination_genome) {
    // NOTE: assumes that both have allocated memory
    // TODO: add error detection for memory allocation?
    int bit_pos;
    if ((*source_genome).genome_len != (*destination_genome).genome_len) {
        fprintf(stderr, "ERROR in SGA_Genome_Copy: Source and Destination of Unequal Length\n");
        exit(0);
    }
 
    (*destination_genome).fitness = (*source_genome).fitness;
    for (bit_pos=0; bit_pos<(*destination_genome).genome_len; bit_pos++) {
        (*destination_genome).bit[bit_pos] = (*source_genome).bit[bit_pos];
    }
}

// Copy source population to dest. population (source unmodified)
void SGA_Population_Copy(population_t *source_population, population_t *destination_population) {
    // TODO: add error detection for genome size within populations
    int p_count;
    if ((*source_population).member_count != (*destination_population).member_count) {
        fprintf(stderr, "ERROR in SGA_Population_Copy: Source and Destination of Unequal Sizes\n");
        exit(0);
    }
    for (p_count=0; p_count < (*source_population).member_count; p_count++) {
        SGA_Genome_Copy(&((*source_population).member[p_count]),&((*destination_population).member[p_count]));
    }
}

void SGA_Population_Metrics(population_t *population, double *metrics) {
    
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
    // (prev. metrics are kept for possible rate convergence calculations - not currently used)
    metrics[5] = metrics[0];
    metrics[6] = metrics[1];
    metrics[7] = metrics[2];
    metrics[8] = metrics[3];
    metrics[9] = metrics[4];

    // find max fitness score and sum of all scores for avg calc
    double max_fitness = population->member[0].fitness;
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

    // decode x, y phenotype of champion
    metrics[3] = SGA_Genome_Decode(&((*population).member[champ_id]), 0, 
                                   ((*population).member[champ_id].genome_len/2)-1, RANGE_MIN, RANGE_MAX);

    metrics[4] = SGA_Genome_Decode(&((*population).member[champ_id]), ((*population).member[champ_id].genome_len/2), 
                                   (*population).member[champ_id].genome_len-1, RANGE_MIN, RANGE_MAX);

    // store avg fitness of generation in metrics
    metrics[1] = (sum_of_fitness / (float)(*population).member_count);

    // Calculate furthest distance between members for diversity metric
    metrics[2] = 0.0;
    for (m_count = 0; m_count < (*population).member_count; m_count++) {
        x1 = SGA_Genome_Decode(&((*population).member[m_count]), 0, 
                               ((*population).member[m_count].genome_len/2)-1, RANGE_MIN, RANGE_MAX);

        y1 = SGA_Genome_Decode(&((*population).member[m_count]), ((*population).member[m_count].genome_len/2), 
                               (*population).member[m_count].genome_len-1, RANGE_MIN, RANGE_MAX);
        
        for (m_inner = m_count+1; m_inner < (*population).member_count; m_inner++) {
            x1 = SGA_Genome_Decode(&((*population).member[m_count]), 0, 
                                   ((*population).member[m_count].genome_len/2)-1, RANGE_MIN, RANGE_MAX);
                               
            y1 = SGA_Genome_Decode(&((*population).member[m_inner]), ((*population).member[m_inner].genome_len/2), 
                                   (*population).member[m_inner].genome_len-1, RANGE_MIN, RANGE_MAX);

            distance = sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0));

            // Store largest distance
            if (distance > metrics[2]) metrics[2] = distance;
        }
    }

    // // calculate % of not unique members for generation (easy but naive algorithm)
    // for (m_count = 0; m_count < (*population).member_count; m_count++) {
    //     for (m_inner = 0; m_inner < m_count; m_inner++) {
    //         // loop through genome bits, exiting if/when a mismatch is found
    //         for (g_bit = 0; g_bit < GENOME_LEN; g_bit++) {
    //             if ((*population).member[m_count].bit[g_bit] != (*population).member[m_inner].bit[g_bit]) {
    //                 break;
    //             }
    //         }
    //         // exit m_inner loop (move onto next iteration of outermost loop) if match found
    //         if (g_bit == GENOME_LEN) break;
    //     }
    //     // increment unique count if loops weren't left early
    //     if (g_bit != GENOME_LEN) num_unique++;
    // }

    // // store non-unique percentage of population
    // metrics[2] = 1.0 - (((float)num_unique) / ((float)POP_SIZE));
}

// ------- Population Initialization Routines -------

// set *already allocated* genome to random bitstring
void SGA_Genome_Init(genome_t *genome, RNG *rng) {
    int bit_count;
    for (bit_count=0; bit_count < (*genome).genome_len; bit_count++) {
        if (rng_uniform01(rng) < 0.5)
            (*genome).bit[bit_count] = 0;
        else
            (*genome).bit[bit_count] = 1;
    }
    (*genome).fitness = 0.0;
}

// randomize all genomes in population
void SGA_Population_Init(population_t *population, RNG *rng) {
    int m_count;
    for (m_count=0; m_count < (*population).member_count; m_count++) {
        SGA_Genome_Init(&((*population).member[m_count]), rng);
    }
}


// ------- Genome Decoding Functions -------

double SGA_Genome_Decode(genome_t *genome, int start, int end, double min, double max) {
    double return_value;
    double max_decimal_value;
    double decimal_value;
    int bit_count;	
    int bit_pos;
  
    decimal_value = 0.0;
    max_decimal_value = 0.0;
    bit_count = 0;

    // Decode the bits from end to start (least to most significant)
    // as a decimal value
    for (bit_pos = end; bit_pos >= start; bit_pos--) {
        decimal_value += (double)((*genome).bit[bit_pos]) * pow(2.0, (double)bit_count);
        max_decimal_value += pow(2.0,(double)bit_count);
        bit_count++;
	}
	
    // calculate and return scaled value between min and max based on decoded decimal
    return min + (max - min) * (decimal_value / max_decimal_value);
}

// void Rosenbrock_Decode(genome_t *genome, double *x, double *y) {
//     const double min_val = -2;
//     const double max_val = 2;
    
//     int start_x = 0;
//     int end_x = (GENOME_LEN / 2) - 1;

//     int start_y = GENOME_LEN / 2;
//     int end_y = GENOME_LEN - 1;

//     *x = SGA_Genome_Decode(genome, start_x, end_x, min_val, max_val);
//     *y = SGA_Genome_Decode(genome, start_y, end_y, min_val, max_val);
// }

// ------- Genome Variation Functions -------

// mutate each bit with probability <mutation_rate>
void SGA_Genome_Mutate(genome_t *genome, double mutation_rate, RNG *rng) {
    // TODO: add error checking for probabilities (between 0 and 1 mutation rate)
    // TODO: make more efficient (replace switch case)
    int bit_count;
    for (bit_count=0; bit_count < (*genome).genome_len; bit_count++) {
        if (rng_uniform01(rng) < mutation_rate) {
            switch ((*genome).bit[bit_count]) {
                case 0: (*genome).bit[bit_count] = 1;
                        break;
                case 1: (*genome).bit[bit_count] = 0;
                        break;
			}
        }
    } 
}

// one point crossover of two genomes
void SGA_Genome_1P_Crossover(genome_t *parent_genome_1, genome_t *parent_genome_2, RNG *rng) {
    // TODO: add error checking for malloc, genome length
    int bit_count;
    int crossover_point;	
    genome_t child_genome_1, child_genome_2;

    // malloc children genomes with same length as parents
    SGA_Genome_Malloc(&child_genome_1, (*parent_genome_1).genome_len);
    SGA_Genome_Malloc(&child_genome_2, (*parent_genome_1).genome_len);
    
    // randomly select crossover_point between 0 and genome_len (uniform distribution)
    crossover_point = (int)round(rng_uniform(rng, 0.0, (double)(*parent_genome_1).genome_len));

    // copy bits before crossover_point to matching child (p1 -> c1, p2 -> c2)
    for (bit_count = 0; bit_count < crossover_point; bit_count++) {
        (child_genome_1).bit[bit_count] = (*parent_genome_1).bit[bit_count];
        (child_genome_2).bit[bit_count] = (*parent_genome_2).bit[bit_count];
    }

    // copy bits including and after crossover_point to opposite child (p1 -> c2, p2 -> c1)
    for (; bit_count < (*parent_genome_1).genome_len; bit_count++) {
        (child_genome_1).bit[bit_count] = (*parent_genome_2).bit[bit_count];
        (child_genome_2).bit[bit_count] = (*parent_genome_1).bit[bit_count];
    }	 
	
    // replace parents with children - GENERATIONAL REPLACEMENT
    SGA_Genome_Copy(&child_genome_1, parent_genome_1);
    SGA_Genome_Copy(&child_genome_2, parent_genome_2); 

    // deallocate temp child memory
    SGA_Genome_Free(&child_genome_1);
    SGA_Genome_Free(&child_genome_2); 
}

// BITWISE MUTATION IMPLEMENTATION
// apply mutation to every genome in population
void SGA_Population_Mutate(population_t *population, double mutation_rate, RNG *rng) {
    // TODO: add error checking?
    int m_count;
	
    // loop through every member of population and apply SGA_Genome_Mutate
    for (m_count = 0; m_count < (*population).member_count; m_count++) {
        SGA_Genome_Mutate(&((*population).member[m_count]), mutation_rate, rng);
    }
}

// ROULETTE WHEEL IMPLEMENTATION
void SGA_Population_Make_New_Generation(population_t *old_population, 
                                        population_t *new_population,
                                        double       mutation_rate,
                                        double       crossover_rate,
                                        double       *metrics,
										RNG          *rng) {

    double sum_of_fitness_scores;
    double *roulette_wheel;
    double max_raw_fitness, min_raw_fitness;
    int m_count;

    int parent_one_index;
    int parent_two_index;
    double die_roll;

    // malloc an array to hold probability density function (PDF) to act as roulette wheel
    roulette_wheel = (double *)malloc(sizeof(double)*(*old_population).member_count); 

    // find max, min fitness scores
    max_raw_fitness = ((*old_population).member[0]).fitness;
    min_raw_fitness = ((*old_population).member[0]).fitness;
	for (m_count = 0; m_count < (*old_population).member_count; m_count++)
	   { if (((*old_population).member[m_count]).fitness > max_raw_fitness)
		     max_raw_fitness = ((*old_population).member[m_count]).fitness;	
		 if (((*old_population).member[m_count]).fitness < min_raw_fitness)
			     min_raw_fitness = ((*old_population).member[m_count]).fitness;
		}	    

    // copy fitness scores into wheel
	for (m_count = 0; m_count < (*old_population).member_count; m_count++)
	    roulette_wheel[m_count] = ((*old_population).member[m_count]).fitness;
	
    // shift scores so that lowest possible fitness is 0
	for (m_count = 0; m_count < (*old_population).member_count; m_count++)
		roulette_wheel[m_count] -= min_raw_fitness;
		
    // ensure not all scores are 0 by adding constant to all
	for (m_count = 0; m_count < (*old_population).member_count; m_count++)
		roulette_wheel[m_count] += 1.0;
		
	// compute sum of all scores, use it to normalize wheel values
	sum_of_fitness_scores = 0.0;
	for (m_count = 0; m_count < (*old_population).member_count; m_count++)
		sum_of_fitness_scores += roulette_wheel[m_count];
	
	for (m_count = 0; m_count < (*old_population).member_count; m_count++)
		roulette_wheel[m_count] /= sum_of_fitness_scores;

    // construct PDF from normalized values
    for (m_count = 1; m_count < (*old_population).member_count; m_count++)
		roulette_wheel[m_count] += roulette_wheel[m_count-1];
		
    // GENERATIONAL REPLACEMENT IMPLEMENTATION LOOP
    for (m_count = 0; m_count < (*old_population).member_count-1; m_count++) {
        // Compute parent indices by looping through PDF until
        // probability >= die_roll is found (once for each parent)
        die_roll = rng_uniform01(rng);
        parent_one_index = 0;
        while (roulette_wheel[parent_one_index] < die_roll) parent_one_index++;
        die_roll = rng_uniform01(rng);
        parent_two_index = 0;
        while (roulette_wheel[parent_two_index] < die_roll) parent_two_index++;
		
        // generate child genomes, mutate, crossover if applicable by rng
        SGA_Genome_Copy( &((*old_population).member[parent_one_index]),&((*new_population).member[m_count]));
        SGA_Genome_Copy( &((*old_population).member[parent_two_index]),&((*new_population).member[m_count+1]));
    
        SGA_Genome_Mutate(&((*new_population).member[m_count]), mutation_rate,rng);
        SGA_Genome_Mutate(&((*new_population).member[m_count+1]), mutation_rate,rng);

        if (rng_uniform01(rng) < crossover_rate)
            SGA_Genome_1P_Crossover(&((*new_population).member[m_count]),&((*new_population).member[m_count+1]),rng);
    }
}


// ------- Genome Fitness Functions -------

// double SGA_Rosenbrock(genome_t *genome) {
//     const double a = 1.0;
//     const double b = 100.0;
    
//     double x;
//     double y;

//     Rosenbrock_Decode(genome, &x, &y);

//     double f_xy = pow((a - x), 2.0) + b * pow((y - x*x), 2.0);

//     // 
//     return 1.0 - f_xy;
// }

double SGA_Himmelblau(genome_t *genome) {
    double x = SGA_Genome_Decode(genome, 0, (genome->genome_len/2)-1, RANGE_MIN, RANGE_MAX);
    double y = SGA_Genome_Decode(genome, (genome->genome_len/2), genome->genome_len-1, RANGE_MIN, RANGE_MAX);

    double f_xy = pow(x*x - 11 + y, 2.0) + pow((x + y*y - 7), 2.0);
    return (0.0 - f_xy);
}


// wrapper for desired fitness function
double SGA_Fitness_Function(genome_t *genome) {
    return SGA_Himmelblau(genome);
}

// apply fitness function wrapped by SGA_Fitness_Function to whole population
void SGA_Population_Compute_Fitness(population_t *population) {
    int m_count;
	
    for (m_count = 0; m_count < (*population).member_count; m_count++)
        ((*population).member[m_count]).fitness = SGA_Fitness_Function(&((*population).member[m_count]));

}

// -----------------------------------------------

int main() {
    
    RNG *random_number_generator;
	population_t POPULATION,POPULATION2;
	int generation_count;
    int converge_count;

    // allocate memory for array to hold running test metrics:
    // max fitness score, avg fitness score, % not unique, champ x, champ y (for generation n and n-1)
    double metrics[10] = {0.0};
	
	// initialize jcg's rng
	random_number_generator = rng_create();
	
	// allocate memory for populations with defined sizes
	SGA_Population_Malloc(&POPULATION, POP_SIZE, GENOME_LEN);  // primary population
	SGA_Population_Malloc(&POPULATION2, POP_SIZE, GENOME_LEN); // "scratch" population 
	
	// initialize/randomize population
	SGA_Population_Init(&POPULATION,random_number_generator);
	
    // SGA_Population_Compute_Fitness(&POPULATION);
    // SGA_Population_Metrics(&POPULATION, metrics);
	// printf("Initial Population\n");
	// SGA_Population_Print(&POPULATION);

	// generation evolutionary loop
	for (generation_count = 0; generation_count < GENERATION_COUNT; generation_count++) {
        // Evaluate Everyone in the population
        SGA_Population_Compute_Fitness(&POPULATION);

        // calculate and print metrics - fitness score of champ, avg fitness score, % identical genomes
        SGA_Population_Metrics(&POPULATION, metrics);
        printf("Himmelblau SGA %d %d %f %f %d %d %f %f %f\n", POP_SIZE, POP_SIZE, MUTATION_RATE, CROSSOVER_RATE,
                                                              generation_count, generation_count*POP_SIZE, 
                                                              metrics[0], metrics[1], metrics[2]);

        //check for convergence - near-optimal solution
        if (metrics[0] > -0.001) break;
        
        // Make the new generation
        SGA_Population_Make_New_Generation(&POPULATION, &POPULATION2, MUTATION_RATE, CROSSOVER_RATE, metrics, random_number_generator);
        // Copy the new population back into the primary population
        SGA_Population_Copy(&POPULATION2, &POPULATION);
    }

    if (generation_count == GENERATION_COUNT) {
        // calculate and print metrics for final generation
        SGA_Population_Compute_Fitness(&POPULATION);
        SGA_Population_Metrics(&POPULATION, metrics);
        printf("Himmelblau SGA %d %d %f %f %d %d %f %f %f\n", POP_SIZE, POP_SIZE, MUTATION_RATE, CROSSOVER_RATE,
                                                              generation_count, generation_count*POP_SIZE, 
                                                              metrics[0], metrics[1], metrics[2]);
    }
		
    printf("\n");	
	
	// printf("Final Population\n");
    // SGA_Population_Print(&POPULATION); 
	
    // if (generation_count < GENERATION_COUNT)
    //     printf("\n[Terminated due to convergence]");
    // else
    //     printf("\n[Termination due to max generation count reached (%d)]", GENERATION_COUNT);
    
    return 0;
}
