#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "CGP-Library-V2.4/src/cgp.h"
#include "tests/cartpole.h"

#define MAX_DEN 10
#define MAX_NEU 100
#define NUM_INPUTS 4
#define NUM_OUTPUTS 1
#define TYPE_OUTPUT 1.0
#define TYPE_NON_OUTPUT -1.0

#define POP_SIZE 20
#define NUM_GEN 50

int neuron_Inputs = 10;
int neuron_Nodes = 15;
int neuron_Outputs = 4;
int neuron_Arity = 2;

int dendrite_Inputs = 10;
int dendrite_Nodes = 15;
int dendrite_Outputs = 4;
int dendrite_Arity = 2;

struct neuron *cur_neuron = NULL;

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

static void copyFunctionSet(struct functionSet *funcSetDest, struct functionSet *funcSetSrc) {

	int i;

	funcSetDest->numFunctions = funcSetSrc->numFunctions;

	for (i = 0; i < funcSetDest->numFunctions; i++) {
		strncpy(funcSetDest->functionNames[i], funcSetSrc->functionNames[i], FUNCTIONNAMELENGTH);
		funcSetDest->functions[i] = funcSetSrc->functions[i];
		funcSetDest->maxNumInputs[i] = funcSetSrc->maxNumInputs[i];
	}
}


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

int nearest_neuron(double pos_x, double pos_y){
    double min_dis = 10;
    int num = 0;
    for(int i = 0; i < MAX_NEU; i++){
        if(cur_neuron[i].valid && cur_neuron[i].pos_x < pos_x){
            if(sqrt(pow((pos_x - cur_neuron[i].pos_x), 2)+pow((pos_y - cur_neuron[i].pos_y), 2)) < min_dis){
                num = i;
            }
        }
    }
    return num;
}

int update_soma(struct neuron n, struct chromosome *neuron_chromo, double performance, double theta_nb, double theta_nd){
    if(!n.valid){
        return -1;
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
    n.bias = outputs[1];
    n.pos_x = outputs[2];
    n.pos_y = outputs[3];
    if(n.health > theta_nb){
        return 1;
    }else if(n.health < theta_nd){
        n.valid = false;
        return -1;
    }else{
        return 0;
    }
}

void update_dendrites(struct neuron n, struct chromosome *dendrite_chromo, double performance, double theta_db, double theta_dd){
    srand(getpid()*time(NULL));
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
            outputs[i] = getChromosomeOutput(dendrite_chromo, j);
        }
        n.d_health[i] = outputs[0];
        n.d_weight[i] = outputs[1];
        n.d_pos_x[i] = outputs[2];
        n.d_pos_y[i] = outputs[3];
        if(n.d_health[i] > theta_db){
            for(int j = 0; j < MAX_DEN; j++){
                if(!n.d_valid[j]){
                    n.d_health[j] = outputs[0];
                    n.d_weight[j] = outputs[1];
                    n.d_valid[j] = true;
                    n.d_pos_x[j] = outputs[2] + -0.1 + (rand() / (double) RAND_MAX * (n.pos_x - -0.1));
                    n.d_pos_y[j] = outputs[3] + -0.1 + (rand() / (double) RAND_MAX * (n.pos_y - -0.1));
                    break;
                }
            }
        }else if(n.d_health[i] < theta_dd){
            n.d_valid[i] = false;
        }
    }
}

int NDS_whi = 1;
int N_ep = 8;
int Num_Problems = 1;
int Epoch_Learning = 1;
double theta_nd_pre = -0.6;
double theta_nb_pre = 0.2;
double theta_dd_pre = -0.7;
double theta_db_pre = 0.2;
double theta_nd_whi = -0.4;
double theta_nb_whi = 0.2;
double theta_dd_whi = -0.65;
double theta_db_whi = -0.6;

void brain_updates(struct neuron *brain_neurons, struct chromosome *neuron_chromo, struct chromosome *dendrite_chromo, 
    double prev_fitness, double theta_nb, double theta_nd, double theta_db, double theta_dd){
    for (int j = 0; j < MAX_NEU; j++){
                int result = update_soma(brain_neurons[j], neuron_chromo, prev_fitness, theta_nb, theta_nd);
                if(result == 0){
                    update_dendrites(brain_neurons[j], dendrite_chromo, prev_fitness, theta_db, theta_dd);
                }else if(result == 1 && brain_neurons[j].type != 1){
                    struct neuron new_n;
                    new_n.health = brain_neurons[j].health;
                    new_n.bias = brain_neurons[j].bias;
                    new_n.pos_x = brain_neurons[j].pos_x -0.1 + (rand() / (double) RAND_MAX * (0.1 - -0.1));
                    new_n.pos_y = brain_neurons[j].pos_y -0.1 + (rand() / (double) RAND_MAX * (0.1 - -0.1));
                    new_n.valid = true;
                    new_n.type = brain_neurons[j].type;
                    for(int k = 0; k < MAX_DEN; k++){
                        new_n.d_valid[k] = brain_neurons[j].d_valid[k];
                        new_n.d_health[k] = brain_neurons[j].d_health[k];
                        new_n.d_weight[k] = brain_neurons[j].d_weight[k];
                        new_n.d_pos_x[k] = brain_neurons[j].d_pos_x[k] -0.1 + (rand() / (double) RAND_MAX * (0.1 - -0.1));
                        new_n.d_pos_y[k] = brain_neurons[j].d_pos_y[k] -0.1 + (rand() / (double) RAND_MAX * (0.1 - -0.1));
                    }
                    for (int l = 0; l < MAX_NEU; l++){
                        if(!brain_neurons[l].valid){
                            brain_neurons[l] = new_n;
                        }
                    }
                    update_dendrites(brain_neurons[j], dendrite_chromo, prev_fitness, theta_db, theta_dd);
                }
            }
}

double epoch_fitness(struct chromosome *neuron_chromo, struct chromosome *dendrite_chromo){
    srand(getpid()*time(NULL));
    struct neuron brain_neurons[MAX_NEU];
    for (int i = 0; i < 100; i++){
        brain_neurons[i].valid = false;
        for (int j = 0; j < MAX_DEN; j++){
            brain_neurons[i].d_valid[j] = false;
        }
    }
    cur_neuron = brain_neurons;

    int NDS_pre = 6;

    //printf("brain initialisation...\n");

    for (int i = 0; i < NUM_OUTPUTS; i++){
        //populate initial output neurons
        brain_neurons[i].bias = 0;
        brain_neurons[i].health = 0.1;
        brain_neurons[i].pos_x = -1 + (rand() / (double) RAND_MAX * (1 - -1));
        brain_neurons[i].pos_y = -1 + (rand() / (double) RAND_MAX * (1 - -1));
        brain_neurons[i].valid = true;
        brain_neurons[i].type = TYPE_OUTPUT;
        for (int j = 0; j < NUM_INPUTS; j++){
            brain_neurons[i].d_health[j] = 0.1;
            brain_neurons[i].d_pos_x[j] = -1 + (rand() / (double) RAND_MAX * (brain_neurons[i].pos_x - -1));
            brain_neurons[i].d_pos_y[j] = -1 + (rand() / (double) RAND_MAX * (brain_neurons[i].pos_y - -1));
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
            brain_updates(brain_neurons, neuron_chromo, dendrite_chromo, prev_fitness, theta_nb_whi, theta_nd_whi, 
                theta_db_whi, theta_dd_whi);
        }

        //printf("extract ANN No. %d\n", ep);
        double fitness = 0;
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
                    /* allocate memory for the node's inputs and connection weights */
                    n->inputs = (int*)malloc(MAX_DEN * sizeof(int));
                    n->weights = (double*)malloc(MAX_DEN * sizeof(double));
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
                    for (int k = 0; k < actarity; k++){
                        n->inputs[k] = inputs[k];
                        n->weights[k] = weights[k];
                    }
                    n->actArity = actarity;
                    n->maxArity = actarity;
                    n->output = 0;
                    extracted->nodes[nodecounter] = n;
                    nodecounter++;
                    if(brain_neurons[j].type == 1){
                        extracted->outputNodes[outputcounter] = nodecounter;
                        outputcounter++;
                    }
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

            printChromosome(extracted, 1);
            struct cartpole_info info;
            reset(&info);

            for (int j = 0; j < 100; j++){
                executeChromosome(extracted, info.state);
                int action = 0;
                if(getChromosomeOutput(extracted, 0) > 1){
                    action = 1;
                }
                double ep_fit = step(&info, action);
                fitness += ep_fit;
            }

            freeChromosome(extracted);
            freeParameters(params);
        }
        prev_fitness = fitness;
        double epoch_fitness = prev_fitness;
        if(Epoch_Learning){
            //update fitness
            if(epoch_fitness < prev_epoch_fitness){
                break;
            }else{
                prev_epoch_fitness = epoch_fitness;
            }
        }
    }
    return prev_epoch_fitness;
}

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
    brain_evaluation();
    return 0;
}
