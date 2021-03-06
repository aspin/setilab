#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <sched.h>
#include <unistd.h>
#include <pthread.h>
#include <math.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

int num_threads;
int num_procs;
pthread_t *tid;

typedef struct {
	int band;
	double bandwidth;
	int band_block;
	int filter_order;
	double* filter_coeffs;
	signal* sig;
	signal* output;
	double* band_power;
} thread_data;

void usage() 
{
    printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands\n");
}

double avg_power(double *data, int num)
{
    int i;
    double ss;
    
    ss=0;
    for (i=0;i<num;i++) { 
	ss += data[i]*data[i];
    }
    
    return ss/num;
}

double max_of(double *data, int num)
{
    double m=data[0];
    int i;
    
    for (i=1;i<num;i++) { 
	if (data[i]>m) { m=data[i]; } 
    }
    return m;
}

double avg_of(double *data, int num)
{
    double s=0;
    int i;
    
    for (i=0;i<num;i++) { 
	s+=data[i];
    }
    return s/num;
}

void remove_dc(double *data, int num)
{
  int i;
  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (i=0;i<num;i++) {
    data[i] -= dc;
  }
}

void *worker(void *arg) {

	// Parse arguments
	thread_data *argcast = (thread_data*) arg;

	// Make the filter
	int i;
	for(i=(argcast->band); i==(argcast->band)+(argcast->band_block); i++) {}
	generate_band_pass(argcast->sig->Fs, 
			   (argcast->band)*((argcast->bandwidth)+0.0001), // keep within limits
			   (argcast->band+1)*((argcast->bandwidth)-0.0001),
			   (argcast->filter_order), 
			   &(argcast->filter_coeffs[(argcast->band)*((argcast->filter_order)+1)]));
	hamming_window((argcast->filter_order),&(argcast->filter_coeffs[(argcast->band)*((argcast->filter_order)+1)]));

	// Convolve
	convolve(argcast->sig->num_samples,
		 argcast->sig->data,
		 argcast->filter_order,
		 &(argcast->filter_coeffs[(argcast->band)*((argcast->filter_order)+1)]),
		 argcast->output->data);

	// Capture characteristics
	(argcast->band_power[argcast->band]) = avg_power(argcast->output->data, argcast->output->num_samples);

	pthread_exit(NULL);

}

int analyze_signal(signal *sig, int filter_order, int num_bands, double *lb, double *ub, int num_procs, int num_threads)
{
    double Fc, bandwidth;
    double filter_coeffs[(filter_order+1)*num_bands];
    double signal_power;
    double band_power[num_bands];
	signal *output[num_bands];
	long tc;
	long rc;

    double start, end;
    
    unsigned long long tstart, tend;
    
    resources rstart, rend, rdiff;
    
    tid = (pthread_t *) malloc(sizeof(pthread_t)*num_threads);

    if (!tid) {
    fprintf(stderr, "cannot allocate memory\n");
    exit(-1);
    }

    Fc=(sig->Fs)/2;
    bandwidth = Fc / num_bands;

    int i;
    for(i=0; i<num_bands; i++) {
    	output[i] = allocate_signal(sig->num_samples, sig->Fs, 0);
	    if (!output[i]) { 
		printf("Out of memory\n");
		return 0;
	    }
    }

    remove_dc(sig->data,sig->num_samples);

    signal_power = avg_power(sig->data,sig->num_samples);

    printf("signal average power:     %lf\n", signal_power);

    get_resources(&rstart,THIS_PROCESS);
    start=get_seconds();
    tstart = get_cycle_count();

    int thread_count = (num_threads < num_bands) ? num_threads : num_bands;
    int band_block = (num_threads < num_bands) ? ceil(num_bands/num_threads): 1;
    int band;
	for (band=0;band<thread_count;band++) {
	if (band == thread_count-1) {
		band_block = num_bands - (band_block * band);
	}
	thread_data datas = {band, bandwidth, band_block, filter_order, &(filter_coeffs[(filter_order+1)*band*band_block]), sig, output[band], band_power};
	rc=pthread_create( &(tid[band]), // thread id gets put here
		       NULL,      // use default attributes
		       worker,    // thread will begin in this function
		       &datas  // we'll give it i as the argument
	             );
	if (rc!=0) { 
	  perror("Failed to start thread");
	  exit(-1);
	}
	}

	// now we will join all the threads
	for (band=0;band<num_threads;band++) {
	rc=pthread_join(tid[band],NULL);   // 
	if (rc!=0) { 
	perror("join failed");
	exit(-1);
	}
	}


    tend = get_cycle_count();
    end = get_seconds();
    get_resources(&rend,THIS_PROCESS);

    get_resources_diff(&rstart, &rend, &rdiff);

    // Pretty print results
    double max_band_power = max_of(band_power,num_bands);
    double avg_band_power = avg_of(band_power,num_bands);
    int wow=0;

#define MAXWIDTH 40

#define THRESHOLD 2.0

#define ALIENS_LOW   50000.0
#define ALIENS_HIGH  150000.0

    *lb=*ub=-1;

    for (band=0;band<num_bands;band++) { 
      double band_low = band*bandwidth+0.0001;
      double band_high = (band+1)*bandwidth-0.0001;
      
      printf("%5d %20lf to %20lf Hz: %20lf ", 
	     band, band_low, band_high, band_power[band]);
      
      for (i=0;i<MAXWIDTH*(band_power[band]/max_band_power);i++) {
	printf("*");
      }
      
      if ( (band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
	   (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) { 

	// band of interest

	if (band_power[band] > THRESHOLD * avg_band_power) { 
	  printf("(WOW)");
	  wow=1;
	  if (*lb<0) { *lb=band*bandwidth+0.0001; }
	  *ub = (band+1)*bandwidth-0.0001;
	} else {
	  printf("(meh)");
	}
      } else {
	printf("(meh)");
      }
      
      printf("\n");
    }

    printf("Resource usages:\n"
	   "User time        %lf seconds\n"
	   "System time      %lf seconds\n"
	   "Page faults      %ld\n"
	   "Page swaps       %ld\n"
	   "Blocks of I/O    %ld\n"
	   "Signals caught   %ld\n"
	   "Context switches %ld\n",
	   rdiff.usertime,
	   rdiff.systime,
	   rdiff.pagefaults,
	   rdiff.pageswaps,
	   rdiff.ioblocks,
	   rdiff.sigs,
	   rdiff.contextswitches);
	   

    printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\nNote that cycle count only makes sense if the thread stayed on one core\n", tend-tstart, cycles_to_seconds(tend-tstart), timing_overhead());
    printf("Analysis took %lf seconds by basic timing\n", end-start);

    return wow;
	
}

int main(int argc, char *argv[])
{
    signal *sig;
    double Fs;
    char sig_type;
    char *sig_file;
    int filter_order;
    int num_bands;
    double start, end;
    
    if (argc!=6) { 
	usage();
	return -1;
    }
    
    sig_type = toupper(argv[1][0]);
    sig_file = argv[2];
    Fs = atof(argv[3]);
    filter_order = atoi(argv[4]);
    num_bands = atoi(argv[5]);
    num_procs = atoi(argv[6]);
    num_threads = atoi(argv[7]);

    assert(Fs>0.0);
    assert(filter_order>0 && !(filter_order & 0x1));
    assert(num_bands>0);

    printf("type:     %s\n"
	   "file:     %s\n"
	   "Fs:       %lf Hz\n"
	   "order:    %d\n"
	   "bands:    %d\n"
	   "processes: %d\n"
	   "threads: %d\n",
	   sig_type=='T' ? "Text" : sig_type=='B' ? "Binary" : sig_type=='M' ? "Mapped Binary" : "UNKNOWN TYPE",
	   sig_file,
	   Fs,
	   filter_order,
	   num_bands,
	   num_procs,
	   num_threads);
    
    printf("Load or map file\n");
    
    switch (sig_type) {
	case 'T':
	    sig = load_text_format_signal(sig_file);
	    break;

	case 'B':
	    sig = load_binary_format_signal(sig_file);
	    break;

	case 'M':
	    sig = map_binary_format_signal(sig_file);
	    break;
	    
	default:
	    printf("Unknown signal type\n");
	    return -1;
    }
    
	if (!sig) { 
	printf("Unable to load or map file\n");
	return -1;
	}

    sig->Fs=Fs;

    if (analyze_signal(sig,filter_order,num_bands,&start,&end, num_procs, num_threads)) { 
	printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n",start,end,(end+start)/2.0);
    } else {
	printf("no aliens\n");
    }

    free_signal(sig);

    return 0;
}


    
