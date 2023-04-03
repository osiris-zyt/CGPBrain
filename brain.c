#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "CGP-Library-V2.4/src/cgp.h"
#include "tests/cartpole.h"

/* Multitasking ANN model project
Author: Yintong Zhang */
/* For original model and details, see https://direct.mit.edu/isal/proceedings/isal2020/473/98437 */

//Global Parameters
#define MAX_DEN 60
#define MAX_NEU 30
#define NUM_INPUTS 4
#define NUM_OUTPUTS 1
#define TYPE_OUTPUT 1.0
#define TYPE_NON_OUTPUT -1.0
#define NDS_whi 1
#define ND_init 5
#define N_ep 8
#define Num_Problems 1
int Epoch_Learning = 1;

#define POP_SIZE 10
#define NUM_GEN 30
#define MAXRECORDS 1000

int neuron_Inputs = 10;
int neuron_Nodes = 15;
int neuron_Outputs = 4;
int neuron_Arity = 2;

int dendrite_Inputs = 10;
int dendrite_Nodes = 15;
int dendrite_Outputs = 4;
int dendrite_Arity = 2;

int seed_global = 0;

struct neuron *cur_neuron = NULL;

double class_dataset[MAXRECORDS][5];

struct neuron{
    double health;
    double bias;
    double pos_x;
    double pos_y;
    int type;
    bool valid;
    double d_health[MAX_DEN];
    double d_weight[MAX_DEN];
    double d_pos_x[MAX_DEN];
    double d_pos_y[MAX_DEN];
    bool d_valid[MAX_DEN];
};

//Singleton random number generator
int getRandom(){
    srand(seed_global);
    seed_global = rand();
    return rand();
}

/*Helper function from CGP Library */
static void copyFunctionSet(struct functionSet *funcSetDest, struct functionSet *funcSetSrc) {

	int i;

	funcSetDest->numFunctions = funcSetSrc->numFunctions;

	for (i = 0; i < funcSetDest->numFunctions; i++) {
		strncpy(funcSetDest->functionNames[i], funcSetSrc->functionNames[i], FUNCTIONNAMELENGTH);
		funcSetDest->functions[i] = funcSetSrc->functions[i];
		funcSetDest->maxNumInputs[i] = funcSetSrc->maxNumInputs[i];
	}
}

/* Custom chromosome node functions for CGP */
double usr_step(const int numInputs, const double *inputs, const double *connectionWeights){
    if(inputs[0] < 0){
        return 0;
    }else{
        return 1;
    }
}

double usr_istep(const int numInputs, const double *inputs, const double *connectionWeights){
    if(inputs[0] < 0){
        return 1;
    }else{
        return 0;
    }
}


int class_prob(struct chromosome *chromo, double **datasets, double *targets, int param_num, int dataset_length){
    int score = 0;
    for(int i = 0; i < dataset_length; i++){
        for(int j = 0; j < param_num; j++){

        }
    }
}

/* Fuction for finding nearest neuron  from a dendrite in an ANN topological graph */
int nearest_neuron(double pos_x, double pos_y){
    double min_dis = 100;
    int num = 0;
    for(int i = 0; i < MAX_NEU; i++){
        if(cur_neuron[i].valid && cur_neuron[i].pos_x < pos_x){
            double dis = sqrt(pow((pos_x - cur_neuron[i].pos_x), 2)+pow((pos_y - cur_neuron[i].pos_y), 2));
            if(dis < min_dis){
                num = i;
                min_dis = dis;
            }
        }
    }
    if(min_dis == 100){
        num = (int)((pow(pos_y, 2))*100)%NUM_INPUTS;
    }
    return num;
}

/* Soma updating function with the soma program */
struct neuron update_soma(struct neuron n, struct chromosome *neuron_chromo, double performance){
    if(!n.valid){
        return n;
    }
    double sum_health, sum_weight, sum_pos_x, sum_pos_y, count;
    for(int i = 0; i < MAX_DEN; i++){
        if(n.d_valid[i]){
            sum_health += n.d_health[i];
            sum_weight += n.d_weight[i];
            sum_pos_x += n.d_pos_x[i];
            sum_pos_y += n.d_pos_y[i];
            count ++;
        }
    }

    double inputs[] = {n.type, n.health, n.bias, n.pos_x, n.pos_y, sum_health/count, sum_weight/count, sum_pos_x/count, 
    sum_pos_y/count, performance};
    executeChromosome(neuron_chromo, inputs);
    double outputs[neuron_Outputs];
    for(int i = 0; i < neuron_Outputs; i++){
        outputs[i] = getChromosomeOutput(neuron_chromo, i);
    }
    n.health = outputs[0];
    n.bias = tanh(outputs[1]);
    n.pos_x = tanh(outputs[2]);
    n.pos_y = tanh(outputs[3]);
    return n;
}

/* Dendrite updating function with the dendrite program */
struct neuron update_dendrites(struct neuron n, struct chromosome *dendrite_chromo, double performance, double theta_db, double theta_dd){
    for(int i = 0; i < MAX_DEN; i++){
        if(!n.d_valid[i]){
            continue;
        }
        int near_num = nearest_neuron(n.d_pos_x[i], n.d_pos_y[i]);
        double inputs[] = {n.bias, n.pos_x, n.pos_y, n.d_health[i], n.d_weight[i], n.d_pos_x[i], n.d_pos_y[i], 
        cur_neuron[near_num].pos_x, cur_neuron[near_num].pos_y, performance};
        executeChromosome(dendrite_chromo, inputs);
        double outputs[dendrite_Outputs];
        for(int j = 0; j < dendrite_Outputs; j++){
            outputs[j] = getChromosomeOutput(dendrite_chromo, j);
        }
        n.d_health[i] = outputs[0];
        n.d_weight[i] = outputs[1];
        n.d_pos_x[i] = tanh(outputs[2]);
        n.d_pos_y[i] = tanh(outputs[3]);
        if(n.d_health[i] > theta_db){
            for(int j = 0; j < MAX_DEN; j++){
                if(!n.d_valid[j]){
                    n.d_health[j] = 0.1;
                    n.d_weight[j] = outputs[1];
                    n.d_valid[j] = true;
                    n.d_pos_x[j] = outputs[2] + -0.1 + (getRandom() / (double) RAND_MAX * (n.pos_x - -0.1));
                    n.d_pos_y[j] = outputs[3] + -0.1 + (getRandom() / (double) RAND_MAX * (n.pos_y - -0.1));
                    break;
                }
            }
        }else if(n.d_health[i] < theta_dd){
            n.d_valid[i] = false;
        }
    }
    return n;
}

//predefined training parameters
double theta_nd_pre = -0.6;
double theta_nb_pre = 0.2;
double theta_dd_pre = -0.7;
double theta_db_pre = 0.2;
double theta_nd_whi = -0.4;
double theta_nb_whi = 0.2;
double theta_dd_whi = -0.65;
double theta_db_whi = -0.6;

/* Brain updating function for a single developmental step */
void brain_updates(struct neuron *brain_neurons, struct chromosome *neuron_chromo, struct chromosome *dendrite_chromo, 
    double prev_fitness, double theta_nb, double theta_nd, double theta_db, double theta_dd){
    for (int j = 0; j < MAX_NEU; j++){
        if(brain_neurons[j].valid){
            brain_neurons[j] = update_soma(brain_neurons[j], neuron_chromo, prev_fitness);
            if(brain_neurons[j].health > theta_nd){
                brain_neurons[j] = update_dendrites(brain_neurons[j], dendrite_chromo, prev_fitness, theta_db, theta_dd);
            }
        }
    }
    for (int j = 0; j < MAX_NEU; j++){
        if(brain_neurons[j].health > theta_nb){
            brain_neurons[j].health = 0.1;

            struct neuron new_n;
            new_n.health = 0.1;
            new_n.bias = brain_neurons[j].bias;
            new_n.pos_x = brain_neurons[j].pos_x + 0.0001;
            new_n.pos_y = brain_neurons[j].pos_y + 0.0001;
            new_n.valid = true;
            new_n.type = TYPE_NON_OUTPUT;
            for(int k = 0; k < MAX_DEN; k++){
                new_n.d_valid[k] = brain_neurons[j].d_valid[k];
                new_n.d_health[k] = brain_neurons[j].d_health[k];
                new_n.d_weight[k] = brain_neurons[j].d_weight[k];
                new_n.d_pos_x[k] = -1 + (getRandom() / (double) RAND_MAX * (new_n.pos_x - -1));
                new_n.d_pos_y[k] = -1 + (getRandom() / (double) RAND_MAX * (1 - -1));
            }
            for (int l = 0; l < MAX_NEU; l++){
                if(!brain_neurons[l].valid){
                    brain_neurons[l] = new_n;
                    break;
                }
            }
        }else if(brain_neurons[j].health < theta_nd && brain_neurons[j].type != TYPE_OUTPUT){
            brain_neurons[j].valid = false;
        }
    }
}

/* ANN development and fitness evaluation for a single epoch */
double epoch_fitness(struct chromosome *neuron_chromo, struct chromosome *dendrite_chromo){
    struct neuron brain_neurons[MAX_NEU];
    for (int i = 0; i < MAX_NEU; i++){
        brain_neurons[i].valid = false;
        for (int j = 0; j < MAX_DEN; j++){
            brain_neurons[i].d_valid[j] = false;
        }
    }
    cur_neuron = brain_neurons;

    int NDS_pre = 6;

    for (int i = 0; i < NUM_OUTPUTS; i++){
        //populate initial output neurons
        brain_neurons[i].bias = 0;
        brain_neurons[i].health = 0.1;
        brain_neurons[i].pos_x = -1 + (getRandom() / (double) RAND_MAX * (1 - -1));
        brain_neurons[i].pos_y = -1 + (getRandom() / (double) RAND_MAX * (1 - -1));
        brain_neurons[i].valid = true;
        brain_neurons[i].type = TYPE_OUTPUT;
        for (int j = 0; j < ND_init; j++){
            brain_neurons[i].d_health[j] = 0.1;
            brain_neurons[i].d_pos_x[j] = -1 + (getRandom() / (double) RAND_MAX * (brain_neurons[i].pos_x - -1));
            brain_neurons[i].d_pos_y[j] = -1 + (getRandom() / (double) RAND_MAX * (1 - -1));
            brain_neurons[i].d_weight[j] = 1.0;
            brain_neurons[i].d_valid[j] = true;
        }
    }

    //printf("pre brain updates...\n");
    for (int i = 0; i < NDS_pre; i++){
        //soma and dendrite update
        brain_updates(brain_neurons, neuron_chromo, dendrite_chromo, 0, theta_nb_pre, theta_nd_pre, 
            theta_db_pre, theta_dd_pre);
    }
    double prev_epoch_fitness = 0;

    //printf("whi brain updates...\n");
    for (int ep = 0; ep < N_ep; ep++){
        double prev_fitness = 0;
        for (int i = 0; i < NDS_whi; i++){
            //evolved soma and dendrite update
            brain_updates(brain_neurons, neuron_chromo, dendrite_chromo, prev_epoch_fitness/100, theta_nb_whi, theta_nd_whi, 
                theta_db_whi, theta_dd_whi);
        }

        //printf("extract ANN No. %d\n", ep);
        // for (int j = 0; j < 25; j++){
        //     printf("%d ", brain_neurons[j].valid);
        // }
        // printf("\n");
        
        double fitness = 0, fitness2 = 0;
        for (int i = 0; i < Num_Problems; i++){
            //extract ANN
            struct chromosome *extracted;
            struct parameters *params;
            params = initialiseParameters(NUM_INPUTS, MAX_NEU, NUM_OUTPUTS, MAX_DEN);
            addNodeFunction(params, "tanh");
            extracted = (struct chromosome*)malloc(sizeof(struct chromosome));
            extracted->nodes = (struct node**)malloc(MAX_NEU * sizeof(struct node*));
            extracted->outputNodes = (int*)malloc(NUM_OUTPUTS * sizeof(int));
            extracted->activeNodes = (int*)malloc(MAX_NEU * sizeof(int));
	        extracted->outputValues = (double*)malloc(NUM_OUTPUTS * sizeof(double));

            int outputcounter = 0, nodecounter = 0;
            for (int j = 0; j < MAX_NEU; j++){
                if(brain_neurons[j].valid){
                    struct node *n;
                    n = (struct node*)malloc(sizeof(struct node));
                    /* set the node's function */
                    n->function = 0;
                    /* set as active by default */
                    n->active = 1;
                    int actarity = 0;
                    double weights[MAX_DEN], inputs[MAX_DEN];
                    for (int k = 0; k < MAX_DEN; k++){
                        if(brain_neurons[j].d_valid[k]){
                            int den_connect = nearest_neuron(brain_neurons[j].d_pos_x[k], brain_neurons[j].d_pos_y[k]);
                            inputs[actarity] = den_connect;
                            weights[actarity] = brain_neurons[j].d_weight[k];
                            actarity++;
                        }
                    }
                    /* allocate memory for the node's inputs and connection weights */
                    n->inputs = (int*)malloc(MAX_DEN * sizeof(int));
                    n->weights = (double*)malloc(MAX_DEN * sizeof(double));
                    for (int k = 0; k < actarity; k++){
                        n->inputs[k] = inputs[k];
                        n->weights[k] = weights[k];
                    }
                    n->actArity = actarity;
                    n->maxArity = actarity;
                    n->output = 0;
                    extracted->nodes[nodecounter] = n;
                    //extracted->activeNodes[nodecounter] = nodecounter+NUM_INPUTS;
                    if(brain_neurons[j].type == 1){
                        extracted->outputNodes[outputcounter] = nodecounter+NUM_INPUTS;
                        outputcounter++;
                    }
                    nodecounter++;
                }
            }
            extracted->arity = MAX_DEN;
            extracted->numNodes = nodecounter;
            extracted->numActiveNodes = nodecounter;
            extracted->numOutputs = outputcounter;
            extracted->numInputs = NUM_INPUTS;
            extracted->nodeInputsHold = (double*)malloc(MAX_DEN * sizeof(double));
            extracted->funcSet = (struct functionSet*)malloc(sizeof(struct functionSet));
	        copyFunctionSet(extracted->funcSet, params->funcSet);

            //removeInactiveNodes(extracted);
            //printChromosome(extracted, 0);
            //setActiveNodes(extracted);

            //Perform a single RL task: CartPole
            struct cartpole_info info;
            reset(&info);
            
            for (int j = 0; j < 1000; j++){
                double ANN_inputs[NUM_INPUTS];
                for(int k = 0; k < 4; k++){
                    ANN_inputs[k] = info.state[k];
                }
                // for(int k = 0; k < 4; k++){
                //     ANN_inputs[k+4] = class_dataset[j][k];
                // }
                executeChromosome(extracted, ANN_inputs);
                int action = 0;
                double target = 0;
                if(getChromosomeOutput(extracted, 0) > 0){
                    action = 1;
                }
                // if(getChromosomeOutput(extracted, 1) > 0){
                //     target = 1;
                // }
                double ep_fit = step(&info, action);
                fitness += ep_fit;
                // if(target == class_dataset[j][4]){
                //     fitness2 += 1;
                // }
            }
            
            freeChromosome(extracted);
            freeParameters(params);
        }
        prev_fitness = (fitness+fitness2);
        double epoch_fitness = prev_fitness;
        if(Epoch_Learning){
            //update fitness
            if(epoch_fitness < prev_epoch_fitness){
                break;
            }else{
                prev_epoch_fitness = epoch_fitness;
            }
        }
        //printf("Fitness for problem 1 is: %f, problem 2 is: %f\n", fitness, fitness2);
    }
    return prev_epoch_fitness;
}

/* Initialize the chromosome population and Evolutionary Algorithm */
void brain_evaluation(){
    struct parameters *neuron_params = NULL;
    struct parameters *dendrite_params = NULL;

    struct chromosome *neuron_chromo[POP_SIZE];
    struct chromosome *dendrite_chromo[POP_SIZE];
    
    neuron_params = initialiseParameters(neuron_Inputs, neuron_Nodes, neuron_Outputs, neuron_Arity);
    dendrite_params = initialiseParameters(dendrite_Inputs, dendrite_Nodes, dendrite_Outputs, dendrite_Arity);

    addNodeFunction(neuron_params, "add,mul,xor");
    addCustomNodeFunction(neuron_params, usr_step, "ustep", -1);
    addCustomNodeFunction(neuron_params, usr_istep, "istep", -1);
    addNodeFunction(dendrite_params, "add,mul,xor");
    addCustomNodeFunction(dendrite_params, usr_step, "ustep", -1);
    addCustomNodeFunction(dendrite_params, usr_istep, "istep", -1);

    for (int i = 0; i < POP_SIZE; i++){
        neuron_chromo[i] = initialiseChromosome(neuron_params);
        dendrite_chromo[i] = initialiseChromosome(dendrite_params);
    }

    printf("start evaluating fitnesses... \n\n");

    double max_fitness = 0;
    int max_number = 0;
    for (int num = 0; num < NUM_GEN; num++){
        for (int i = 0; i < POP_SIZE; i++){
            double ep_fitness = epoch_fitness(neuron_chromo[i], dendrite_chromo[i]);
            printf("The epoch fitness is: %f\n", ep_fitness);
            if(ep_fitness > max_fitness){
                max_fitness = ep_fitness;
                max_number = i;
            }
        }
        printf("The max fitness in generation %d is: %f\n", num+1, max_fitness);
        for (int i = 0; i < POP_SIZE; i++){
            copyChromosome(neuron_chromo[i], neuron_chromo[max_number]);
            copyChromosome(dendrite_chromo[i], dendrite_chromo[max_number]);
        }
        for (int i = 0; i < POP_SIZE; i++){
            mutateChromosome(neuron_params, neuron_chromo[i]);
            mutateChromosome(dendrite_params, dendrite_chromo[i]);
        }
    }


    //free everything
    for (int i = 0; i < POP_SIZE; i++){
        freeChromosome(neuron_chromo[i]);
        freeChromosome(dendrite_chromo[i]);
    }
    freeParameters(neuron_params);
    freeParameters(dendrite_params);
    return 0;
}

int main(void) {
    seed_global = 114514;

    //Reading data for a classification problem (Used only in solving more than 1 problems)
    FILE *file;
    int nRecords = 0;
    memset(class_dataset, 0, sizeof(class_dataset));
    file = fopen("banknote_auth.txt", "r");

    if (file == NULL){
        printf("File does not exist!");
    }else{
        //printf("before\n");
        while (EOF != fscanf(file, "%f,%f,%f,%f,%f", &class_dataset[nRecords][0], &class_dataset[nRecords][1], 
                                &class_dataset[nRecords][2], &class_dataset[nRecords][3], &class_dataset[nRecords][4]) && nRecords<MAXRECORDS)
        {
            //printf("%d\n", nRecords);
            nRecords++;
        }
    }

    fclose(file);

    brain_evaluation();

    /* struct cartpole_info info;
    reset(&info);
    double ep_fit;
    ep_fit += step(&info, 0);
    printf("%f\n", info.state[2]);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 0);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 0);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 0);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 0);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    printf("%f\n", info.state[2]);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
    ep_fit += step(&info, 1);
    printf("%f\n", ep_fit);
 */
    return 0;
}
