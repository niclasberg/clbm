#ifndef PLOT_H_
#define PLOT_H_

typedef struct {
	double min_val;
	double max_val;
} PlotOptions;

void destroy_plot();
void init_plot();
void imagesc(unsigned int, unsigned int, double *, PlotOptions *);


#endif /* PLOT_H_ */
