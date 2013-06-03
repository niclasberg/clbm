#ifndef WORKERPOOL_H_
#define WORKERPOOL_H_
#include <string.h>


typedef void (* job_predicate_t)(void *);

void workerpool_run();
void workerpool_init(size_t);
void workerpool_destroy();
void workerpool_push_job(job_predicate_t, void *);


#endif /* WORKERPOOL_H_ */
