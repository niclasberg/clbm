#include "workerpool.h"
#include <stdio.h>
#include <time.h>

void test_func(void * arg) {
	int id = *((int *) arg);

	printf("In thead %d\n", id);
	sleep(1);
	printf("Thead %d done\n", id);
}

int main(int argc, char ** argv) {
	int i;
	int thread_data[100];

	workerpool_init(10);

	for(i=0; i < 100; ++i) {
		thread_data[i] = i;
		workerpool_push_job(test_func, (void *) &thread_data[i]);
	}

	workerpool_run();
}
