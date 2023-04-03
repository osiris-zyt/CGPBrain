#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include "cartpole.h"

/* A C version of the classic CartPole Problem. 
Original code (in Python): https://github.com/openai/gym/blob/master/gym/envs/classic_control/cartpole.py 
Author: Yintong Zhang */

#define GRAVITY 9.8
#define MASSCART 1
#define MASSPOLE 0.1
#define TOTAL_MASS 1.1
#define LENGTH 0.5
#define POLEMASS_LENGTH 0.05
#define FORCE_MAG 10
#define TAU 0.02

double theta_threshold_radians = 12*2*M_PI/360;
double x_threshold = 2.4;
int steps_beyond_terminated = 0;


double * returnState(struct cartpole_info info){
    return info.state;
}

void reset (struct cartpole_info *info) {
    srand(steps_beyond_terminated);
    double min = -0.5;
    double max = 0.5;
    double min_rad = -0.1;
    double max_rad = 0.1;
    for(int i = 0; i < 2; i++){
        info->state[i] = min + (rand() / (double) RAND_MAX * (max - min));
    }
    for(int i = 2; i < 4; i++){
        info->state[i] = min_rad + (rand() / (double) RAND_MAX * (max_rad - min_rad));
    }
    info->reward = 0;
    info->terminated = false;
    steps_beyond_terminated = info->state[0]*time(NULL);
}

double step (struct cartpole_info *info, int action){
    double x = info->state[0];
    double x_dot = info->state[1];
    double theta = info->state[2];
    double theta_dot = info->state[3];
    double costheta = cos(theta);
    double sintheta = sin(theta);

    double force = 0;
    if(action == 1){
        force = FORCE_MAG;
    }else {
        force = -FORCE_MAG;
    }
    
    double temp = (force + POLEMASS_LENGTH * pow(theta_dot, 2) * sintheta) / TOTAL_MASS;
    double thetaacc = (GRAVITY * sintheta - costheta * temp) / (LENGTH *
     (4.0 / 3.0 - MASSPOLE * pow(costheta, 2) / TOTAL_MASS));
    double xacc = temp - POLEMASS_LENGTH * thetaacc * costheta / TOTAL_MASS;

    x = x + TAU * x_dot;
    x_dot = x_dot + TAU * xacc;
    theta = theta + TAU * theta_dot;
    theta_dot = theta_dot + TAU * thetaacc;

    info->state[0] = x;
    info->state[1] = x_dot;
    info->state[2] = theta;
    info->state[3] = theta_dot;

    bool terminated = info->terminated;
    //printf("threshold: %f, %f\n", x_threshold, theta_threshold_radians);
    if(x < -x_threshold || x > x_threshold || theta < -theta_threshold_radians || theta > theta_threshold_radians){
        terminated = true;
    }

    double reward = 0.0;
    if(!terminated){
        reward = 1.0;
    // }else if (steps_beyond_terminated == NULL){
    //     steps_beyond_terminated = 0;
    //     reward = 1.0;
    }else {
        //if(steps_beyond_terminated == 0){
            //printf("calling 'step()' even when current environment is terminated.\n");
        //}
        //steps_beyond_terminated += 1;
        reward = 0.0;
    }

    info->reward = reward;
    info->terminated = terminated;

    return reward;
}
