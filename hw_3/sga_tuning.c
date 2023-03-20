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
		
// function prototypes to keep GCC happy with the ordering of definitions and calls
double SGA_Genome_Decode(genome_t *, int, int, double, double);
double SGA_Gray_Decode(genome_t *, int, int, double, double);

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
    for (int member_c=0; member_c < (*population).member_count; member_c++) {
	    SGA_Genome_Free(&((*population).member[member_c]));
    }
    
    free((*population).member);
    (*population).member_count = 0;
}


// ------- Population Utility Routines -------

// print out genome as a bitstring with current fitness score in parenthesis
void SGA_Genome_Print(genome_t *genome, double range_min, double range_max){
    int bit_count;
    double x = SGA_Gray_Decode(genome, 0, (genome->genome_len/2)-1, range_min, range_max);
    double y = SGA_Gray_Decode(genome, (genome->genome_len/2), genome->genome_len-1, range_min, range_max);

    for (bit_count=0; bit_count < (*genome).genome_len; bit_count++) {
        printf("%d",(*genome).bit[bit_count]);
    }
    printf(" [%f, %f] (%lf)\n", x, y, (*genome).fitness);
	
}

// Print ALL the genomes in a population, with fitness scores
void SGA_Population_Print(population_t *population, double range_min, double range_max) {
    int p_count;
    
    for (p_count=0; p_count < (*population).member_count; p_count++) {
        SGA_Genome_Print(&((*population).member[p_count]), range_min, range_max);
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

void SGA_Population_Metrics(population_t *population, double *metrics, double range_min, double range_max) {
    
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
    metrics[3] = SGA_Gray_Decode(&((*population).member[champ_id]), 0, 
                                   ((*population).member[champ_id].genome_len/2)-1, range_min, range_max);

    metrics[4] = SGA_Gray_Decode(&((*population).member[champ_id]), ((*population).member[champ_id].genome_len/2), 
                                   (*population).member[champ_id].genome_len-1, range_min, range_max);

    // store avg fitness of generation in metrics
    metrics[1] = (sum_of_fitness / (float)(*population).member_count);

    // Calculate furthest distance between members for diversity metric
    metrics[2] = 0.0;
    // DISABLED FOR PERFORMANCE
    // for (m_count = 0; m_count < (*population).member_count; m_count++) {
    //     x1 = SGA_Gray_Decode(&((*population).member[m_count]), 0, 
    //                            ((*population).member[m_count].genome_len/2)-1, range_min, range_max);

    //     y1 = SGA_Gray_Decode(&((*population).member[m_count]), ((*population).member[m_count].genome_len/2), 
    //                            (*population).member[m_count].genome_len-1, range_min, range_max);
        
    //     for (m_inner = m_count+1; m_inner < (*population).member_count; m_inner++) {
    //         x1 = SGA_Gray_Decode(&((*population).member[m_count]), 0, 
    //                                ((*population).member[m_count].genome_len/2)-1, range_min, range_max);
                               
    //         y1 = SGA_Gray_Decode(&((*population).member[m_inner]), ((*population).member[m_inner].genome_len/2), 
    //                                (*population).member[m_inner].genome_len-1, range_min, range_max);

    //         distance = sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0));

    //         // Store largest distance
    //         if (distance > metrics[2]) metrics[2] = distance;
    //     }
    // }

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

// Treat bitstring as a gray code - convert to equivalent binary then decode as normal
double SGA_Gray_Decode(genome_t *genome, int start, int end, double min, double max) {
    int bit;	
    int bit_pos;
    double return_value;

    genome_t genome_temp;

    // Allocate memory for and make a copy of the genome
    SGA_Genome_Malloc(&genome_temp, (*genome).genome_len);
    SGA_Genome_Copy(genome, &genome_temp);

    // printf("\n");
    // for (bit_pos=0; bit_pos < genome_temp.genome_len; bit_pos++) {
    //     printf("%d", genome_temp.bit[bit_pos]);
    // }

    // Convert copy's gray bits to binary from most to least significant
    // b_n      = g_n
    // b_n-1    = g_n XOR g_n-1             = b_n XOR g_n-1
    // b_n-2    = g_n XOR g_n-1 XOR g_n-2   = b_n-1 XOR g_n-2
    // (etc.)
    for (bit_pos = start+1; bit_pos <= end; bit_pos++) {
        // Calculate bit as XOR of current and prev. bit
        // NOTE: casting ASCII chars 1 and 0 gives ints 48 and 49, but XOR still works as expected/desired
        bit = genome_temp.bit[bit_pos] ^ genome_temp.bit[bit_pos-1];
        // printf("%d ^ %d = %d\n", genome_temp.bit[bit_pos], genome_temp.bit[bit_pos-1], bit);
        
        // Convert int 0 or 1 to corresponding ASCII char
        genome_temp.bit[bit_pos] = bit;
	}
	
    // printf("\n");
    // for (bit_pos=0; bit_pos < genome_temp.genome_len; bit_pos++) {
    //     printf("%d", genome_temp.bit[bit_pos]);
    // }

    // Decode the copied and converted genome as a binary bitstring
    return_value = SGA_Genome_Decode(&genome_temp, start, end, min, max);
    SGA_Genome_Free(&genome_temp);
    return return_value;
}


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

// two point crossover of two genomes
void SGA_Genome_2P_Crossover(genome_t *parent_genome_1, genome_t *parent_genome_2, RNG *rng) {
    // TODO: add error checking for malloc, genome length
    int bit_count;
    int crossover_p1, crossover_p2, temp;	
    genome_t child_genome_1, child_genome_2;

    // malloc children genomes with same length as parents
    SGA_Genome_Malloc(&child_genome_1, (*parent_genome_1).genome_len);
    SGA_Genome_Malloc(&child_genome_2, (*parent_genome_1).genome_len);
    
    // randomly select crossover_p1 between 0 and genome_len-1 (uniform distribution)
    crossover_p1 = (int)floor(rng_uniform(rng, 0.0, (double)(*parent_genome_1).genome_len));
    crossover_p2 = crossover_p1;

    // Select necessarily different p2
    while (crossover_p2 == crossover_p1) crossover_p2 = (int)floor(rng_uniform(rng, 0.0, (double)(*parent_genome_1).genome_len));

    // Swap points if necessary so that p1 is always less than p2
    if (crossover_p2 < crossover_p1) {
        temp = crossover_p2;
        crossover_p2 = crossover_p1;
        crossover_p1 = temp;
    }

    // copy bits before crossover_p1 to matching child (p1 -> c1, p2 -> c2)
    for (bit_count = 0; bit_count < crossover_p1; bit_count++) {
        (child_genome_1).bit[bit_count] = (*parent_genome_1).bit[bit_count];
        (child_genome_2).bit[bit_count] = (*parent_genome_2).bit[bit_count];
    }

    // copy bits between p1 and p2, including p1, to opposite child (p1 -> c2, p2 -> c1)
    for (; bit_count < crossover_p2; bit_count++) {
        (child_genome_1).bit[bit_count] = (*parent_genome_2).bit[bit_count];
        (child_genome_2).bit[bit_count] = (*parent_genome_1).bit[bit_count];
    }	 

    // copy bits after p2 to matching child (p1 -> c1, p2 -> c2)
    for (; bit_count < (*parent_genome_1).genome_len; bit_count++) {
        (child_genome_1).bit[bit_count] = (*parent_genome_1).bit[bit_count];
        (child_genome_2).bit[bit_count] = (*parent_genome_2).bit[bit_count];
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
                                        int          selection_pwr,
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
	for (m_count = 0; m_count < (*old_population).member_count; m_count++) {
	    if (((*old_population).member[m_count]).fitness > max_raw_fitness) {
            max_raw_fitness = ((*old_population).member[m_count]).fitness;
        }
		if (((*old_population).member[m_count]).fitness < min_raw_fitness)
            min_raw_fitness = ((*old_population).member[m_count]).fitness;
    }	    

    // copy fitness scores into wheel
	for (m_count = 0; m_count < (*old_population).member_count; m_count++)
	    roulette_wheel[m_count] = ((*old_population).member[m_count]).fitness;
	
    // Scale scores to be in (0, 1] based on min and max
	for (m_count = 0; m_count < (*old_population).member_count; m_count++) {
		roulette_wheel[m_count] = (roulette_wheel[m_count] - min_raw_fitness) / (max_raw_fitness - min_raw_fitness) + 0.001;
        continue;
    }
		
	// compute sum of all scores, use it to normalize wheel values
	sum_of_fitness_scores = 0.0;
	for (m_count = 0; m_count < (*old_population).member_count; m_count++) {
        // SELECTION PRESSURE EXPERIMENT
        for (int i = 0; i < (selection_pwr - 1); i++) {
            roulette_wheel[m_count] *= roulette_wheel[m_count];
        }
        sum_of_fitness_scores += roulette_wheel[m_count];
    }
	
	for (m_count = 0; m_count < (*old_population).member_count; m_count++)
		roulette_wheel[m_count] /= sum_of_fitness_scores;

    // construct CDF from normalized values
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

    free(roulette_wheel);
}


// ------- Genome Fitness Functions -------

double SGA_Rosenbrock(genome_t *genome, double range_min, double range_max) {
    const double a = 1.0;
    const double b = 100.0;
    
    double x = SGA_Gray_Decode(genome, 0, (genome->genome_len/2)-1, range_min, range_max);
    double y = SGA_Gray_Decode(genome, (genome->genome_len/2), genome->genome_len-1, range_min, range_max);

    double f_xy = pow((a - x), 2.0) + b * pow((y - x*x), 2.0);

    return 0.0 - f_xy;
}

double SGA_Himmelblau(genome_t *genome, double range_min, double range_max) {
    double x = SGA_Gray_Decode(genome, 0, (genome->genome_len/2)-1, range_min, range_max);
    double y = SGA_Gray_Decode(genome, (genome->genome_len/2), genome->genome_len-1, range_min, range_max);

    double f_xy = pow(x*x - 11 + y, 2.0) + pow((x + y*y - 7), 2.0);
    return (0.0 - f_xy);
}


// wrapper for desired fitness function
double SGA_Fitness_Function(genome_t *genome, double range_min, double range_max) {
    return SGA_Rosenbrock(genome, range_min, range_max);
}

// apply fitness function wrapped by SGA_Fitness_Function to whole population
void SGA_Population_Compute_Fitness(population_t *population, double range_min, double range_max) {
    int m_count;
	
    for (m_count = 0; m_count < (*population).member_count; m_count++)
        ((*population).member[m_count]).fitness = SGA_Fitness_Function(&((*population).member[m_count]), range_min, range_max);

}

// -----------------------------------------------

int main() {
    
    // ------- (Initial/Default) EA Parameters -------

    int pop_size            = 100;
    int genome_len          = 20;

    double mutation_rate    = 0.005;
    double mutation_ctrl    = 0.9;
    double crossover_rate   = 0.925;
    int selection_pwr       = 2;

    double range_min        = -5.12;
    double range_max        = 5.11;


    // ------- (TODO) Test Parameters -------

    int generation_count    = 500;
    double optimum_fitness  = -0.001;
    int run_count           = 100;

    double mutation_min     = 0.4;
    double mutation_max     = 0.6;
    double mutation_ctrl_min     = 0.00001;
    double mutation_ctrl_max     = 0.009;
    double crossover_min    = 0.5;
    double crossover_max    = 1.0;
    double selection_max    = 5.0;
    int step_count          = 5;        // Number of values tested for each parameter INCLUDING end points

    // -------


    RNG *random_number_generator;
	population_t POPULATION,POPULATION2;
	int generation;
    double mutation_mod;            // Used for time-based mutation control
    double generation_avg = 0.0;    // Stores running sum of generations to reach optimum (or not) for avg calculation
    double champ_avg = 0.0;         // ^^^ Same for best fitness
    double champ_stddev = 0.0;
    FILE *csv;

    // allocate memory for array to hold running test metrics:
    // max fitness score, avg fitness score, % not unique, champ x, champ y (for generation n and n-1)
    double metrics[10] = {0.0};
	
	// initialize jcg's rng
	random_number_generator = rng_create();
	
	// allocate memory for populations with defined sizes
	SGA_Population_Malloc(&POPULATION, pop_size, genome_len);  // primary population
	SGA_Population_Malloc(&POPULATION2, pop_size, genome_len); // "scratch" population 

    for (mutation_rate = mutation_min; mutation_rate < (mutation_max*1.0001); mutation_rate += ((mutation_max - mutation_min)/(double)(step_count-1))) {
        printf("[mutation_rate = %f]\n", mutation_rate);

        for (mutation_ctrl = mutation_ctrl_min; mutation_ctrl < (mutation_ctrl_max*1.0001); mutation_ctrl += ((mutation_ctrl_max - mutation_ctrl_min)/(double)(step_count-1))) {
            printf("\t[mutation_ctrl = %f]\n", mutation_ctrl);

            for (crossover_rate = crossover_min; crossover_rate < (crossover_max*1.0001); crossover_rate += ((crossover_max - crossover_min)/(double)(step_count-1))) {
                printf("\t\t[crossover_rate = %f]\n", crossover_rate);

                for (selection_max = 1.0; selection_max < 9.0001; selection_max += ((9.0 - 1.0)/(double)(step_count-1))) {
                    printf("\t\t\t[selection_max = %f]\n", selection_max);

                    generation_avg = 0.0;
                    champ_stddev = 0.0;
                    champ_avg = 0.0;
                    for (int i = 0; i < run_count; i++) {

                        // initialize/randomize population
                        SGA_Population_Init(&POPULATION,random_number_generator);
                        
                        // SGA_Population_Compute_Fitness(&POPULATION, range_min, range_max);
                        // SGA_Population_Metrics(&POPULATION, metrics, range_min, range_max);
                        // printf("Initial Population\n");
                        // SGA_Population_Print(&POPULATION, range_min, range_max);
                    
                        // generation evolutionary loop
                        for (generation = 0; generation < generation_count; generation++) {
                            // Evaluate Everyone in the population
                            SGA_Population_Compute_Fitness(&POPULATION, range_min, range_max);

                            // calculate and print metrics - fitness score of champ, avg fitness score, % identical genomes
                            SGA_Population_Metrics(&POPULATION, metrics, range_min, range_max);
                            // printf("Rosenbrock SGA %d %d %f %f %d %d %f %f %f\n", pop_size, pop_size, mutation_rate, crossover_rate,
                            //                                                       generation, generation*pop_size, 
                            //                                                       metrics[0], metrics[1], metrics[2]);

                            // check for convergence - near-optimal solution
                            // if (metrics[0] > optimum_fitness) break;
                            
                            // Calculate mutation modifier (reduces mutation each generation) - TIME-BASED CONTROL
                            // mutation_mod decreases linearly, multiplied with mutation_rate to give linear decay over time (generations)
                            // mutation_ctrl defines the minimum value of mutation_mod
                            mutation_mod = ((double)generation_count - (double)generation*(1.0 - mutation_ctrl)) / (double)generation_count;

                            selection_pwr = 1 + floor(selection_max * (double)generation / (double)generation_count);

                            // Make the new generation
                            SGA_Population_Make_New_Generation(&POPULATION, &POPULATION2, mutation_rate*mutation_mod, crossover_rate, selection_pwr, metrics, random_number_generator);
                            // Copy the new population back into the primary population
                            SGA_Population_Copy(&POPULATION2, &POPULATION);
                        }

                        if (generation == generation_count) {
                            // calculate and print metrics for final generation
                            SGA_Population_Compute_Fitness(&POPULATION, range_min, range_max);
                            SGA_Population_Metrics(&POPULATION, metrics, range_min, range_max);
                            // printf("Rosenbrock SGA %d %d %f %f %d %d %f %f %f\n", pop_size, pop_size, mutation_rate, crossover_rate,
                            //                                                       generation, generation*pop_size, 
                            //                                                       metrics[0], metrics[1], metrics[2]);
                        }

                        // printf("Final Population\n");
                        // SGA_Population_Print(&POPULATION, range_min, range_max); 
                        
                        // if (generation_count < GENERATION_COUNT)
                        //     printf("\n[Terminated due to convergence]");
                        // else
                        //     printf("\n[Termination due to max generation count reached (%d)]", GENERATION_COUNT);

                        generation_avg += generation;
                        champ_stddev += pow(metrics[0], 2.0);
                        champ_avg += metrics[0];
                        // printf("\t\t\t\t[%f, %f]: %f; gen avg = %f\n", metrics[3], metrics[4], metrics[0], metrics[1]);
                    }

                    generation_avg /= run_count;
                    champ_avg /= run_count;
                    champ_stddev = sqrt((champ_stddev / (double)run_count) - pow(champ_avg, 2.0));

                    csv = fopen("sga_results.csv", "a");

                    if (csv == NULL) {
                        printf("Error opening sga_results.csv");
                        exit(1);
                    }
                    fflush(csv);
                    fprintf(csv, "%f, %f, %f, %f, %f, %f\n", mutation_ctrl, mutation_rate, crossover_rate, selection_max, champ_avg, champ_stddev);
                    fclose(csv);

                }
            }
        }
    }
		
    printf("\n");	
    
    return 0;
}
