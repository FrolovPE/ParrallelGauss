#include "lib.h"
class args
{
public:

char *name{};
double *a{};
double *b{};
int *mainblocks{};
double *minnorms{};
int n{};
int m{};
int r{};
int s{};
int thr{};
int p{};
int err{};
double time{};

pthread_barrier_t *barrier{};
pthread_mutex_t *mutex{};


args() = default;
args(const args& x) = delete;
args(args&& x) = default;
args& operator=(const args& x) = delete;
args& operator=(args&& x) = default;
~args() = default;

};