#include "workerpool.h"
#include <string.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

// Threads
pthread_t * threads;
size_t num_threads;

// Job struct
struct job_node {
	job_predicate_t predicate;
	void * arg;
	struct job_node * next;
};

// The jobs are stored in a FIFO type queue
pthread_mutex_t job_mutex;
struct job_node * head_job;
struct job_node * tail_job;
size_t remaining_jobs;

/*
 * push_job
 * Pushes a job to the job queue. The jobs are added at the end.
 */
void workerpool_push_job(job_predicate_t pred, void * arg)
{
	// Create the node to add
	struct job_node * new_node = (struct job_node *) malloc(sizeof(struct job_node));
	new_node->next = NULL;
	new_node->predicate = pred;
	new_node->arg = arg;

	// Lock the list so other threads cannot make concurrent changes
	pthread_mutex_lock(&job_mutex);

	// Push to the job list
	if( ! head_job) {
		// No jobs exist in the queue
		head_job = new_node;
		tail_job = new_node;
	} else {
		tail_job->next = new_node;
		tail_job = new_node;
	}

	++remaining_jobs;

	// Unlock the list
	pthread_mutex_unlock(&job_mutex);
}

void * do_work(void * arg)
{
	struct job_node * current_job;

	while(1) {
		pthread_mutex_lock(&job_mutex);

		// Pop job from queue and update current head
		if(head_job) {
			--remaining_jobs;
			current_job = head_job;
			head_job = head_job->next;
		} else {
			// No jobs remaining, exit the loop
			pthread_mutex_unlock(&job_mutex);
			break;
		}

		/*printf("*** Workerpool: job started, %d jobs remaining in queue. ***\n", (int)remaining_jobs);*/

		pthread_mutex_unlock(&job_mutex);

		// Execute the predicate function and free the job node once done
		current_job->predicate(current_job->arg);
		/*printf("*** Workerpool: job completed ***\n");*/
		free(current_job);
	}

	pthread_exit(NULL);
}

void workerpool_run()
{
	/*printf("*** Workerpool: execution of jobs started ***\n");*/

	size_t i;
	int rc;

	// Create thread attributes
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	// Start all the threads
	for(i = 0; i < num_threads; ++i) {
		rc = pthread_create(&threads[i], &attr, do_work, NULL);

		if(rc) {
			printf("ERROR: pthread_create, error code: %d\n", rc);
			exit(-1);
		}
	}

	pthread_attr_destroy(&attr);

	// Wait for all threads to finish
	for(i = 0; i < num_threads; ++i) {
		rc = pthread_join(threads[i], NULL);
	}
}

void workerpool_init(size_t thread_count)
{
	// Initialize job list
	head_job = NULL;
	tail_job = NULL;
	remaining_jobs = 0;
	pthread_mutex_init(&job_mutex, NULL);

	// Allocate space for threads
	num_threads = thread_count;
	threads = (pthread_t *) malloc(num_threads * sizeof(pthread_t));
}

void workerpool_destroy()
{
	// Destroy mutex
	pthread_mutex_destroy(&job_mutex);

	// Free the allocated threads
	free(threads);
	threads = NULL;
	num_threads = 0;

	// Clear the job queue
	struct job_node * next_job;
	while(head_job) {
		next_job = head_job->next;
		free(head_job);
		head_job = next_job;
	}

	head_job = NULL;
	tail_job = NULL;
	remaining_jobs = 0;
}
