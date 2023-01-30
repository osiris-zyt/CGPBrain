#ifndef HEADER_FILE
#define HEADER_FILE

struct cartpole_info
{
    /* data */
    double state[4];
    double reward;
    bool terminated;

};


extern void reset (struct cartpole_info *info);
extern double step (struct cartpole_info *info, int action);
extern double * returnState(struct cartpole_info info);

#endif


