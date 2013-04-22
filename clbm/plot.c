#include "plot.h"
#include "gnuplot_i.h"

gnuplot_ctrl * plot_handle;
char temp_file_name[] = "temp";
char temp_file_name2[] = "temp2";

void imagesc(unsigned int lx, unsigned int ly, double * u, PlotOptions * plot_opts) {
	//Setup gnuplot
	gnuplot_cmd(plot_handle, "set multiplot layout 1, 2");
	gnuplot_cmd(plot_handle, "unset title");
	gnuplot_cmd(plot_handle, "set pm3d map");
	gnuplot_cmd(plot_handle, "set xrange [%d:%d]", 0, lx-1);
	gnuplot_cmd(plot_handle, "set yrange [%d:%d]", 0, ly-1);
	gnuplot_cmd(plot_handle, "set cbrange [%f:%f]", plot_opts->min_val, plot_opts->max_val);
	gnuplot_cmd(plot_handle, "set palette rgbformulae 22,13,10");

	// Create a temporary file for the data
	FILE * temp_file = fopen(temp_file_name, "w");
	FILE * temp_file2 = fopen(temp_file_name2, "w");

	// Write data
	size_t i, j;

	for(j = 0; j < ly; ++j) {
		for(i = 0; i < lx; ++i) {
			fprintf(temp_file, "%.18e\n", u[i*ly + j]);
		}
		fprintf(temp_file, "\n");
		fprintf(temp_file2, "%.18e\t%.18e\n", u[(lx/2)*ly + j], (double) j);
	}

	fclose(temp_file);
	fclose(temp_file2);

	// Feed the file to gnuplot
	gnuplot_cmd(plot_handle, "splot '%s'", temp_file_name);
	gnuplot_cmd(plot_handle, "unset pm3d");
	gnuplot_cmd(plot_handle, "set xrange [%f:%f]", plot_opts->min_val, plot_opts->max_val);
	gnuplot_cmd(plot_handle, "unset cbrange");
	gnuplot_cmd(plot_handle, "plot '%s' using 1:2 with lines lc rgb '#0025ad'", temp_file_name2);
	gnuplot_cmd(plot_handle, "unset multiplot");

}

void init_plot()
{
	plot_handle = gnuplot_init();
}

void destroy_plot()
{
	gnuplot_close(plot_handle);
}
