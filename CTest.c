/*==================================================================================================================================
This file was created for testing/development for the C programming language							                           | 
CREATOR: Rohan Jyoti                                                                                                               |
FOCUS: Audio Programming and generic C code testing                                                                                |
Inspired by The Audio Programming Book by Richard Boulanger and Victor Lazzarini                                                   |
MLA Citation: (Boulanger, Richard), (Lazzarini, Victor). The Audio Programming Book. Cambridge, MA: The MIT Press, 2011. Print.    |
                                                                                                                                   |
==================================================================================================================================*/

//To Run: ./CTest [function] [function args] (./ can be omitted on Mac OS X)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//#include <my_global.h>
//#include <mysql.h>

#include <portsf.h>
#include <wave_generator.h>
#include <time.h>
#include <breakpoints.h>



#ifndef CTest_h
#define CTest_h
//For extendability purposes, move this definition section to its own header file
void function1();
void iscale(int argc, char *argv[]);
void expdecay(int argc, char *argv[]);
void sinetext(int argc, char *argv[]);
int byte_order();
void sf2float(int argc, char *argv[]);
void sig_gen(int argc, char *argv[]);
void osc_gen(int argc, char *argv[]);
void summation_function(int argc, char *argv[]);
void secant_method(int argc, char *argv[]);
float fnx(float x);
void oddSumTest(int argc, char *argv[]);
void pointerTest();
void aQRec();
float m_POW(int a, int b);
#endif

#ifndef M_PI
//if M_PI is not defined in math.h (cross-compiler || cross-platform...)
#define M_PI (3.141592654)
#endif

//The following enum is for use with function: sinetext
//Basically, it assigns numbers with meaningful variable names
enum {ARG_NAME,ARG_NSAMPS,ARG_FREQ,ARG_SR,ARG_NCHANNELS,ARG_NARGS};
//so ARG_NAME = 0; ARG_NSAMPS = 1; ARG_FREQ = 2; ARG_SR = 3; ARG_NARGS = 4;

//The following enum is for use with function: my_sql_1
enum {NAME,CR_REM,HOST,DBUSER,DBPASSWD,DBNAME,INFILE};

//The following enum is for use with function: sf2float
enum {PROGNAME, IN_FILE, OUTFILE, NUMARGS};

//The following enums are for use with function: sig_gen
enum {SIG_GEN_PROGNAME, SIG_GEN_OUTFILE, SIG_GEN_TYPE, SIG_GEN_DUR, SIG_GEN_SRATE, SIG_GEN_AMP, SIG_GEN_FREQ, SIG_GEN_NARGS};
enum {SINE_WAVE = 0, TRIANGLE_WAVE, SQAURE_WAVE, SAW_UP_WAVE, SAW_DOWN_WAVE, NUM_WAVE_TYPES};
///The following definition is for use with function: sig_gen
#define NFRAMES (1024)

//The following enums are for use with function: osc_gen
enum {OSC_GEN_PROGNAME, OSC_GEN_OUTFILE, OSC_GEN_DUR, OSC_GEN_SRATE, OSC_GEN_NUM_CHANS, OSC_GEN_AMP, OSC_GEN_FREQ, OSC_GEN_WAVETYPE, OSC_GEN_NUM_OSCS, OSC_GEN_NARGS};
enum {OSC_GEN_SQUARE = 0, OSC_GEN_TRIANGLE, OSC_GEN_SAW_UP, OSC_GEN_SAW_DOWN, OSC_GEN_PULSE};


int main(int argc, char *argv[])
{
	if(argc==1) 
    {
        printf("You must enter function\n");
        exit(1);     
    }
    char *fn_name;
    fn_name=argv[1];
    if(strcmp(fn_name,"function1")==0) function1(); //basic pointer aritmetic
    else if(strcmp(fn_name,"iscale")==0) iscale(--argc, ++argv); //generate E.T tables for N-notes to the octave (N<=24) ***NOTICE --argc and ++argv
    else if(strcmp(fn_name,"expdecay")==0) expdecay(--argc, ++argv); //generate exponential attack or decay breakpoint data
    else if(strcmp(fn_name,"sinetext")==0) sinetext(--argc, ++argv); //generate sinusoidal curve for use with GNUplot --> showcases use of enum
    else if(strcmp(fn_name,"byte_order")==0) { if(byte_order()==0) printf("System is big-endian\n"); else printf("System is little-endian\n"); } //tells endian type of system
    else if(strcmp(fn_name,"sf2float")==0) sf2float(--argc, ++argv); //convert soundfile to float format
    else if(strcmp(fn_name,"sig_gen")==0) sig_gen(--argc, ++argv); //generate oscillation signal (sine, triangle, square, sawtooth up, sawtooth down)
	else if(strcmp(fn_name,"osc_gen")==0) osc_gen(--argc, ++argv); //Oscillator Bank (Array of oscillators) for additive sythesis for the natural specturm of sound (classic square, trianlgle, saw_up, saw_down)
	else if(strcmp(fn_name,"summation_function")==0) summation_function(--argc, ++argv); // summation for a^n (CS473 Quiz3)
	else if(strcmp(fn_name,"secant_method")==0) secant_method(--argc, ++argv); // secant method implementation for CS357 Quiz3
	else if(strcmp(fn_name,"oddSumTest")==0) oddSumTest(--argc, ++argv); // oddSum program for Rishi
	else if(strcmp(fn_name,"pointerTest")==0) pointerTest();//pointer test for CS241
	else if(strcmp(fn_name,"aQRec")==0) aQRec(); //tests for audioQueueRecorder
    else
    {
        printf("Function not found\n");
        exit(1);
    }
    
    return 0;
}

void function1()
{
    /*pointer arithmetic --> sets each element in buffer array to 0*/
    int num=16;
    double buffer[num];
    double *ptr = buffer; //Set ptr to the start of buffer array --> prt = buffer[0]
    int i;
    for(i=0; i<num; i++)
    {
        *ptr=0.0; //set the value of current address the ptr is pointing to
        ptr++; //compiler knows to do ptr=ptr+sizeof(double)
    }
    
    //verify
    for(i=0; i<num; i++)
    {
        if(i!=num-1) printf("%d,", (int)buffer[i]);
        else printf("%d.\n", (int)buffer[i]);
    }
	
	//This is a test for CS473 SP12 HW4 Q2
	
	int j, k;
	int n = 4;
	
	for(i = 0; i<n; i++)
	{
		for(j = 0; j<i+1; j++)
		{
			for(k=0; k<j-1; k++)
			{
				printf("(%d,%d) (%d,%d)\n", i, j, j, k);
			}
		}
	}
	
	long int value = 0xFF;
	long int flag = 0XFF;
	if( !(value & ~flag) ) printf ("They are equal\n");
	
	int a = 16;
	int b = 4;
	printf("(%d << %d = %d ;; %d*2^%d = %f\n", a, b, a << b, a, b, a*pow(2,b));
	
	a=2;
	b=3;
	for(i=0; i<10; i++)
		printf("m_POW(%d,%d) = %d\n", a, i, (int)m_POW(a,i));
	
	
	int flag2;
	if(flag2) printf("flag2 is true with value: %d\n", flag2);
	
}

void iscale(int argc, char *argv[])
{
    /* usage--> iscale [-m] [-i] N startval [outfile.txt]
        -m      :   sets format of startval as MIDI note
        -i      :   prints the calcualted interval as well as the abs freq
        outfile :   optional text filename for output data
    */
    
    int notes, i, ismidi=0, write_interval=0, err=0;
    double startval, basefreq, ratio;
    FILE* fp;
    double intervals[25];
    
    //verify
    for(i=0; i<argc; i++)
    {
        printf("argv[%d]: %s\n", i, argv[i]);
    }
    
    //Check first arg for flag option: argc at least 2?
    /*The follwing will make use of pointer arithmetic (with argc and argv)
        as well as make use of the following concept: recall that argv is 
        an array of array (pointer to an array of char array); therefore,
        we can use argv[1][0] to go the first character of the second argument
    */
    while(argc>1)
    {
        if(argv[1][0]=='-')
        {
            //a flag has been sensed
            if(argv[1][1]=='m') ismidi=1;
            else if(argv[1][1]=='i') write_interval=1;
            else
            {
                printf("error: unrecongized option %s\n", argv[1]);
                exit(1);
            }
            //pointer arithmetic
            argc--;
            argv++;
        }
        else break; //meaning that we are passed the '-' flag section and into other parameter which will be handled later        
    }
    
    if(argc<3)
    {
        //NOTICE because of previous pointer aritmetic on argc, argc is now effectively its original value minus at most 2 (for -m and -i)
        //recall that a maximum of 5 arguments is possible --> 5-2=3 <-- therefore, argc should be >= 3
        printf("Insufficient arguments\n");
        printf("Usage --> iscale [-m] [-i] N startval [outfile.txt]\n");
        exit(1);
    }
    
    //Now, again due to previous pointer arithmetic on argv, we now expect argv[1] to hold N and argv[2] to hold startval
    notes = atoi(argv[1]);
    if(notes<1 || notes>24)
    {
        printf("N value out of range; N --> 1<=N<=24\n");
        exit(1);
    }
    startval = atof(argv[2]);
    if(ismidi)
    {
        //meaning -m flag was input --> ismidi=1 --> TRUE
        if(startval<0.0 || startval>127.0)
        {
            printf("startval out of range; startval with midi flag --> 0<=startval<=127\n");
            exit(1);
        }
    }
    else
    {
        //meaning no -m flag --> ismidi=0 --> FALSE --> dealing with hz frequencies
        //Since it is a freq, we must have positive numbers
        if(startval<=0.0)
        {
            printf("startval out of range; startval without midi flag means dealing with frequencies --> startval>0\n");
            exit(1);
        }
        
    }
    
    //Check for optional filename
    fp = NULL;
    if(argc==4)
    {
        //meaning filename is provided
        fp = fopen(argv[3], "w"); //r = read-only; w = write; a = append
        if(fp==NULL)
        {
            printf("WARNING: Unable to create file %s\n", argv[3]);
            perror("");
            exit(1);
        }
    }
    
    //=========================All Parameters Ready --> Fill Array and Write To File if Created=========================
    
    //Find basefreq, if val is MIDI
    if(ismidi)
    {
        double c0, c5;
        //find base MIDI note
        ratio = pow(2.0, 1.0/12.0); //2^(1/12)
        c5 = 220.0 * pow(ratio, 3);
        c0 = c5 * pow(0.5, 5);
        basefreq = c0 * pow(ratio, startval);
    }
    else
    {
        //in the case val is not midi but frequencies
        basefreq = startval;
    }
    
    //calculate ratio from notes, and fill the array
    ratio = pow(2.0, 1.0/notes);
    for(i=0; i<=notes; i++)
    {
        intervals[i] = basefreq;
        basefreq+=ratio;
        
        //Finaly, read the array, write to screen, and optimally to file
        double r_i = pow(ratio, i);
        if(write_interval) printf("%d: \t%f\t%f\n", i, r_i, intervals[i]);
        else printf("%d: \t%f\n", i, intervals[i]);
        
        //write to file
        if(fp)
        {
            if(write_interval) err = fprintf(fp, "%d \t%f\t%f\n", i, r_i, intervals[i]);
            else err = fprintf(fp, "%d: \t%f\n", i, intervals[i]);
            if(err<0)
            {
                printf("file write error\n");
                perror("");
            }
        }
        
    }
    
    //close FILE file descriptor
    if(fp) 
    {
        if(err>=0) printf("Output file %s successfully created\n", argv[3]);
        fclose(fp);
    }
    printf("Task Finished...\n");   
}

void expdecay(int argc, char *argv[])
{
    //This function will generate exponential attack or decay brekapoint data to be used with GNUplot
    //It makes use of STDOUT and STDERR to differentiate output streams so that when data is ready to be used with GNUplot, we can use > redirect to txt file
    //without also redirecting any error messages we may or may not display; > redirect operator takes STDOUT stream and writes it to a file
    
    /********************USAGE*******************
     expdecay duration npoints startval endval
    ********************************************/
    
    //Additional Info --> after output to STDOUT on terminal --> we will convert to dB log scale so it is easier to see
    //Audio amplitude values between 0 and 1 are said to be normalized
    //p(dB) = 20.0 log10(x)
    
    //To use GNUplot --> aquaTerm must be installed
        //At command prompt type gnuplot --> plot "expdecay.txt" with lines
        //For dB log scale --> plot "expdecay.txt" using (20.0 * log10($2+0.00001)) with lines <-- the 0.00001 is to normalize 0 values since log(0) is illegal
    
    int i;
    printf("argc count is %d\n", argc);
    for(i=0; i<argc; i++)
    {
        printf("argv[%d]: %s\n", i, argv[i]);
    }
    
    int npoints;
    double startval, endval;
    double dur, step, start, end, thisstep;
    double fac, valrange, offset;
    const double verysmall = 1.0e-4; // ~ -80dB  //
    
    if(argc!=5)
    {
        fprintf(stderr, "Usage: expdecay duration npoints startval endval\n");
        exit(1);
    }
    
    dur = atof(argv[1]);
    npoints = atoi(argv[2]);
    if(dur<=0 || npoints<=0)
    {
        fprintf(stderr, "error: both duration and npoints must be positve values\n");
        exit(1);
    }
    
    step = dur/npoints;
    startval = atof(argv[3]);
    endval = atof(argv[4]);
    valrange = endval - startval;
    if(valrange==0.0)
    {
        fprintf(stderr, "warning: start and end values are the same\n");
    }
    
    //initialize normalized exponential as attack or decay
        //decay
    if(valrange<0)
    {
        //meaning startval is greater than endval
        start = 1.0;
        end = verysmall;
        valrange = -valrange;
        offset = endval;
    }
    else
    {
        //attack
        //meaning valrange >0 --> startval < endval
        start = verysmall;
        end = 1.0;
        offset = startval;
    }
    
    thisstep = 1.0;
    //make normalized curve, scale output to input values, range
    fac = pow(end/start, 1.0/npoints);
    for(i=0; i<npoints; i++)
    {
        fprintf(stdout, "%.4lf\t%.8lf\n", thisstep, offset + (start * valrange));
        start *= fac;
        thisstep += step;
    }
    //print final value
    fprintf(stdout, "%.4lf\t%.8lf\n", thisstep, offset + (start * valrange));
    
    fprintf(stderr, "task finised\n");
    
}

void sinetext(int argc, char *argv[])
{
    //generates sinusoidal samples for use with GNUplot --> showcases enum
    
    /**************************USAGE*************************
     sinetext nsamps freq samprate nchannels > sine.txt
     example: sinetext 50 440 44100  2 > sine.txt
        50 sample, 440 Hz = Conecrt A, 
        44100 Hz = standard CD sample rate
        1 = mono; 2 = stereo <-- multi-channel implementation
     *******************************************************/
    
    //To plot in GNUplot --> plot "sine.txt" with impulses
    //or for DAC (Digital ot Analog converter) format --> plot "sine.txt" with steps
    //For MULTI-CHANNEL support --> plot "sine.txt" using($1) with lines, "sine.txt" using($2) with lines
    
    int i,j,nsamps,nchannels;
    double samp,freq,samprate;
    double twopi = 2.0 * M_PI;
    double angleincr;
    
    if(argc!=ARG_NARGS)
    {
        //meaning that arg count is not equal to 4 as defined by enum
        fprintf(stderr,"Usage: sinetext nsamps freq srate\n");
        exit(1);
    }
    nsamps = atoi(argv[ARG_NSAMPS]);
    freq = atof(argv[ARG_FREQ]);
    samprate = atof(argv[ARG_SR]);
    nchannels = atoi(argv[ARG_NCHANNELS]);
    
    angleincr = twopi*freq/samprate;
    
    for(i=0; i<nsamps; i++)
    {
        samp = sin(angleincr*i);
        //virtual implementation of mult-channel output
        //in reality, the n-columns of data would
        //be meaningful, here we are just using sqaure-wave
        //methodolgy for multi-channel effect
        for(j=1; j<=nchannels; j++)
        {
            if(j==nchannels) fprintf(stdout, "%.8lf\n", pow(samp, j));
            else fprintf(stdout, "%.8lf\t", pow(samp, j)); //sample^(1,2,3,...,nchannels)
        }
        //fprintf(stdout, "%.8lf\n", samp);
    }
    
    fprintf(stderr, "Task Finished\n");
}

int byte_order()
{
    //This is a trick acquired from SNDAN programmers to tell endianness of system
    //0 = big-endian machine; 1 = little-endian machine
    /*
     decimal        :   1164413355
     hexadecimal    :   0x456789AB
     big-endian     :   45  67  89  AB
     little-endian  :   AB  89  67  45
     
     each column is a byte; Intel uses little-endian; Big-endian synonymous with left to right reading as in English
     
     */
    
    int one = 1;
    char *endptr = (char *)&one;
    return (*endptr);
}



void sf2float(int argc, char *argv[])
{
    //This function will sound portsf to convert a soundfile to float format
    
    /*=========================USAGE=========================
     sf2float infile outfile
     ======================================================*/
    if(argc!=NUMARGS)
    {
        fprintf(stderr, "USAGE: sf2float infile outfile\n");
        exit(1);
    }
    
    
    PSF_PROPS props; //structure containing sample rate, channels, sample type, format, channel format...
    long framesread,totalread;
    int infd =-1 , outfd = -1;
    psf_format outformat = PSF_FMT_UNKNOWN; //initially make the outformat unknown until properly defined by create function
    PSF_CHPEAK* peaks = NULL; //in waveform and aiff, peaks holds the channel peak values accordingly
    float* frame = NULL; //a frame is a set of n datatypes (16-bit --> short; 24 bit --> long; 32 bit --> float) where n is determined by number of channels
    
    //initialize portsf
    int psfd = psf_init();
    if(psfd<0)
    {
        fprintf(stderr, "Unable to initialize portsf\n");
        exit(1);
    }
    
    //function to open exsiting soudfile protoype--> psf_sndOpen(const char *path, PSF_PROPS *props, int rescale)
        //PSF_PROPS <-- refer to note above on definition line of props
        //rescale --> normally set to 0 for samples to be read unaltered; otherwise, rescale normalizes sound between -1 and 1 to preven clipping
    infd = psf_sndOpen(argv[IN_FILE], &props, 0);
    if(infd<0)
    {
        fprintf(stderr, "Unable to open soundfile\n");
        goto exit; //Works like a JAL command; exit block must be within same function block; used for cleanup purposes. ALL SUBSEQUENT ERROR HANDLING MUST USE GOTO
    }
    
    //tell user if in_file is already floats (i.e 32-bit)
    if(props.samptype == PSF_SAMP_IEEE_FLOAT)
    {
        fprintf(stderr, "Infile already of type float\t Nothing to do...\n");
        goto exit;
    }
    else props.samptype = PSF_SAMP_IEEE_FLOAT; //if not already of type float --> then set the samptype parameter of props strcuture to type float
    
    //Check to see that outfile extension is applicable --> noncompressed, lossless --> wav, aiff, aif, afc, aifc
    outformat = psf_getFormatExt(argv[OUTFILE]);
    if(outformat==PSF_FMT_UNKNOWN)
    {
        fprintf(stderr, "Improper outfile format\nAccepted Formats: .wav, .aiff, .aif, .afc, .aifc\n");
        goto exit; 
    }
    else props.format=outformat; //if format applicable --> set the format paramter of props structure to the outformat as discovered by portfs on outfile extension provided
    
    //Opening a sound file for writing prototype --> psf_sndCreate(const char *path, const PSF_PROPS *props, int clip_floats, int minheader, int mode)
        //clip_floats --> set the way the floats are written <-- used in conjunction with rescale value being nonzero in psf_sndOpen
        //minheader --> normally set to 0; set to 1 for legacy compatibility
        //mode --> control of read-write access to the file {PSF_CREATE_RDWR, PSF_CREATE_TEMPORARY, PSF_CREATE_WRONLY}
    outfd = psf_sndCreate(argv[OUTFILE], &props, 0, 0, PSF_CREATE_RDWR);
    if(outfd<0)
    {
        fprintf(stderr, "Unable to create soundfile\n");
        goto exit;
    }
    
    //Allocate space for one sample frame
    frame = (float*)malloc(props.chans * sizeof(float));
    if(frame==NULL) goto no_memory;
    //Allocate space for PEAK info
    peaks = (PSF_CHPEAK*)malloc(props.chans * sizeof(PSF_CHPEAK));
    if(peaks==NULL) goto no_memory;
    
    //single-frame loop --> copying as floats (conversion process)
    framesread = psf_sndReadFloatFrames(infd, frame, 1);
    totalread = 0;
    while(framesread==1)
    {
        totalread++;
        if(psf_sndWriteFloatFrames(outfd,frame,1) != 1)
        {
            fprintf(stderr, "Error writing to file\n");
            break;
        }
        /* ===============DO ANY PROCESSING HERE; CORE OF AUDIO PROCESSING===============*/
        framesread = psf_sndReadFloatFrames(infd,frame,1);
    }
    if(framesread<0)
    {
        fprintf(stderr, "Error handling infile; outfile is incomplete\n");
        goto exit;
    }
    else fprintf(stdout, "Task Finished %d sample frames copied to %s\n", (int)totalread, argv[OUTFILE]);
    
    //report PEAK values to user
    if(psf_sndReadPeaks(outfd, peaks, NULL) >0)
    {
        long i;
        double peakTime;
        fprintf(stdout, "PEAK Information:\n");
        for(i=0; i<props.chans; i++)
        {
            peakTime = (double) peaks[i].pos / props.srate;
            fprintf(stdout, "CH %d:\t%.4f at %.4f secs\n", (int)(i+1), peaks[i].val, peakTime);
        }
    }
    
    exit:
    if(infd>=0) psf_sndClose(infd);
    if(outfd>=0) psf_sndClose(outfd);
    if(frame) free(frame);
    if(peaks) free(peaks);
    psf_finish(); // close portfs
    fprintf(stderr, "Program Terminated\n");
    exit(1);
    
    no_memory:
    fprintf(stderr, "No Memory!\n");
    goto exit;
}

void sig_gen(int argc, char *argv[])
{
    //This function will generate simple waveforms --> Sine, Triangle, Square, Sawtooth up, Sawtooth down
	//Showcases Pseudo-Object Oriented Programming, Function Pointers (functoids)
    /***************USAGE***************
     siggen outfile wavetype dur srate amp freq
     **********************************/
    if(argc < SIG_GEN_NARGS)
    {
        fprintf(stderr, "insufficient arguments.\n"
               "usage: siggen outfile wavetype dur srate amp freq\n"
               "where wavetype =:\n"
               "                0 = sine\n"
               "                1 = triangle\n"
               "                2 = square\n"
               "                3 = sawtooth up\n"
               "                4 = sawtooth down\n"
               "dur   = duration of outfile (seconds)\n"
               "srate = required sample rate of outfile\n"
               "amp   = amplitude value or breakpoint file (0 < amp <= 1.0)\n"
               "freq  = frequency value (freq > 0) or breakpoint file.\n"
               );
        exit(1);
    }
    
   /********* initialize all dynamic resources to default states **********/
    PSF_PROPS outprops;
    int ofd = -1;
    PSF_CHPEAK* peaks = NULL;
    psf_format outformat = PSF_FMT_UNKNOWN; //all psf* defined in portsf sound library
    long nframes = NFRAMES;
    float* outframe = NULL;
    double amp,freq,dur;
    unsigned long outframes,nbufs,i;
    long remainder;
    int wavetype;
    OSCIL* osc = NULL; //OSCIL pseudo-object defined in wave_generator.h
    node_func_ptr node; //function pointer defined in wave_generator.h
    FILE* fpfreq = NULL;
    FILE* fpamp = NULL;
    BRKSTREAM* freqstream = NULL; //BRKSTREAM* defined in breakpoint.h
    BRKSTREAM* ampstream = NULL;
    unsigned long brkfreqSize = 0, brkampSize = 0;
    double minval, maxval;
    
    /********** error checking on wave type **********/
    wavetype = atoi(argv[SIG_GEN_TYPE]);
    if(wavetype < SINE_WAVE || wavetype >= NUM_WAVE_TYPES)
    {
        fprintf(stderr, "Error: wavetype not defined\n");
        goto exit;
    }
    
	/********** assigning function pointer accordingly **********/
    switch(wavetype)
    {
        case(SINE_WAVE):
            node = sine_wave_node; //node is a function pointer (type: node_func_ptr); sine_wave_node is function
            break;
		case(SQAURE_WAVE):
			node = square_wave_node;
			break;
		case(TRIANGLE_WAVE):
			node = triangle_wave_node;
			break;
		case(SAW_UP_WAVE):
			node = sawtooth_up_wave_node;
			break;
		case(SAW_DOWN_WAVE):
			node = sawtooth_down_wave_node;
			break;
    }
    
    /********** define outfile format (with embedded error checking) **********/
    outprops.srate = atoi(argv[SIG_GEN_SRATE]);
    if(outprops.srate <= 0)
    {
        fprintf(stderr, "Error: srate must by positive (>=0)\n");
        goto exit;
    }
    outprops.chans = 1; //as of this sig_gen implementation --> mono
    outprops.samptype = (psf_stype) PSF_SAMP_16; //16-bit samples
    outprops.chformat = STDWAVE; //STDWAVE defined in portsf library 
    
    dur = atof(argv[SIG_GEN_DUR]);
    if(dur <= 0.0)
    {
        fprintf(stderr, "Error: duration must be positive (>=0)\n");
        goto exit;
    }
    outframes = (unsigned long) (dur * outprops.srate + 0.5);
    nbufs = outframes / nframes;
    remainder = outframes - nbufs * nframes;
    if(remainder) nbufs++;
    
    /********** open breakpoint files, or set constants if no breakpoint file is defined **********/
    /***** amp (amplitude) breakpoints file *****/
    fpamp = fopen(argv[SIG_GEN_AMP], "r"); //try opening specified argument
    if(fpamp == NULL)
    {
        //meaning that either argument SIG_GEN_AMP is not a breakpoint file or error opening
        amp = atof(argv[SIG_GEN_AMP]);
        if(amp <= 0.0 || amp > 1)
        {
            fprintf(stderr, "Error: amplitude must be between 0 and 1 (0 <= amp <= 1)\n");
            goto exit;
        }
    }
    else
    {
        //meaning that argument SIG_GEN_AMP is indeed a breakpoint file
        ampstream = bps_newstream(fpamp,outprops.srate,&brkampSize);
        if(ampstream == NULL)
        {
            fprintf(stderr, "Error reading amplitude breakpoint file: %s\n", argv[SIG_GEN_AMP]);
            goto exit;
        }
        if(bps_getminmax(ampstream, &minval, &maxval))
        {
            //meaning we called upon said function that will put into minval and maxval appropriate values
            //if said function returns 1 --> we have error
            fprintf(stderr, "Error reading range breakpoint file: %s\n", argv[SIG_GEN_AMP]);
            goto exit;
        }
        if(minval < 0.0 || minval > 1.0 || maxval < 0.0 || maxval > 1.0)
        {
            fprintf(stderr, "Error: amplitude values out of range in breakpoint file\n");
            goto exit;
        }
    }
    
    /***** freq (frequency) breakpoints file *****/
    fpfreq = fopen(argv[SIG_GEN_FREQ], "r");
    if(fpfreq == NULL)
    {
        //meaning that either argument SIG_GEN_FREQ is not a breakpoint file or error opening
        freq = atof(argv[SIG_GEN_FREQ]);
        if(freq <= 0.0)
        {
            fprintf(stderr, "Error: frequency must be positive\n");
            goto exit;
        }
    }
    else
    {
        //meaning that argument SIG_GEN_FREQ is indeed a breakpoint file
        freqstream = bps_newstream(fpfreq, outprops.srate, &brkfreqSize);
        if(freqstream == NULL)
        {
            fprintf(stderr, "Error reading frequncy from breakpoint file: %s\n", argv[SIG_GEN_FREQ]);
            goto exit;
        }
        if(bps_getminmax(freqstream, &minval, &maxval))
        {
            fprintf(stderr, "Error reading range from breakpoint file: %s\n", argv[SIG_GEN_FREQ]);
            goto exit;
        }
        if(minval <= 0.0)
        {
            fprintf(stderr, "Error in breakpoing file; Frequency values must be positive\n");
            goto exit;
        }
            
    }
    
    
    /********** Initialize OSCIL pseudo-object **********/
    osc = new_oscil(outprops.srate);
    if(osc == NULL)
    {
        fprintf(stderr, "No memory for oscillator\n");
        goto exit;
    }
    
    /********** Initialize portsf **********/
    if(psf_init())
    {
        fprintf(stderr, "Error: Unable to initialize portsf\n");
        goto exit;
    }
    
    /********** Sample Frame Buffer Memory Aloocation **********/
    outframe = (float*)malloc(nframes * outprops.chans * sizeof(float));
    if(outframe == NULL)
    {
        fprintf(stderr, "Error: No memory for outframe buffer\n");
        goto exit;
    }
    
    /********** Check for correct output extension **********/
    outformat = psf_getFormatExt(argv[SIG_GEN_OUTFILE]);
    if(outformat == PSF_FMT_UNKNOWN)
    {
        fprintf(stderr, "Error: invalid extension\n"
                "Possible extensions:\t *.wav\t *.aiff\t *.aif\t *.afc\t *.aifc\n");
        goto exit;
    }
    //reach here meaning outformat is viable extension
    outprops.format = outformat;
    
    /********** Initialize Peak nfo **********/
    peaks = (PSF_CHPEAK*)malloc(outprops.chans * sizeof(PSF_CHPEAK));
    if(peaks == NULL)
    {
        fprintf(stderr, "Error: No memory for peaks\n");
        goto exit;
    }
    
    /********** Create and Handle outfile **********/
    ofd = psf_sndCreate(argv[SIG_GEN_OUTFILE],&outprops, 0, 0, PSF_CREATE_RDWR); //refer to sf2float for documentation
    if(ofd < 0)
    {
        fprintf(stderr, "Error: Unable to create outfile\n");
        goto exit;
    }
    
    /********** Audio Processing **********/
    for(i=0; i<nbufs; i++)
    {
        long j;
        if(i == nbufs-1) 
            nframes = remainder;
        
        for(j=0; j<nframes; j++)
        {
            if(freqstream)
                freq = bps_tick(freqstream);
            if(ampstream)
                amp = bps_tick(ampstream);
            outframe[j] = (float)(amp * node(osc,freq)); //call upon function via function pointer with OSCIL pseudo-object and freq;
        }
        if(psf_sndWriteFloatFrames(ofd, outframe, nframes) != nframes)
        {
            fprintf(stderr, "Error: unable to write to outfile\n");
            goto exit;
        }
    }   
    
    /********** Clean_up section **********/
    exit:
    fprintf(stderr, "Task Finished\n");
    /***** Report PEAK values to user *****/
    if(psf_sndReadPeaks(ofd,peaks,NULL) > 0){
		long i;
		double peaktime;
		printf("PEAK information:\n");	
        for(i=0;i < outprops.chans;i++){
            peaktime = (double) peaks[i].pos / (double) outprops.srate;
            printf("CH %d:\t%.4f at %.4f secs\n", (int)(i+1), peaks[i].val, peaktime);
        }
	}
    
    
    if(ofd >= 0)
		if(psf_sndClose(ofd))
			printf("Error closing outfile %s\n",argv[SIG_GEN_OUTFILE]);
	if(outframe)
		free(outframe);
	if(peaks)
		free(peaks);
	/*TODO: cleanup any other resources */
	if(osc)
		free(osc);
	if(fpfreq)
		if(fclose(fpfreq))
			printf("Error closing breakpoint file %s\n",argv[SIG_GEN_FREQ]);
	if(fpamp)
		if(fclose(fpamp))
			printf("Error closing breakpoint file %s\n",argv[SIG_GEN_AMP]);
	if(freqstream){
		bps_freepoints(freqstream);
		free(freqstream);
	}
	if(ampstream){
		bps_freepoints(ampstream);
		free(ampstream);
	}
	psf_finish();
    exit(1);
}

void osc_gen(int argc, char *argv[])
{
	//This function serves as an Osccilator bank (an arrary of oscillators) for additive synthesis
	//It produces more natural sounding waves confined in the spectrum of the availabke ranges of the sampling rate at or above the Nyquist limit
	
	//Structurally, this program is very similar to sig_gen, except that we will only be calling upon sine_wave_node and implement additive synthesis
	//to achieve the other waveforms. If osc_gen_num_oscs == 1 (meaning only 1 oscillator) a normal sine wave will be produced.
	
	//This function makes use of wave_generator.c/.h objects and functions
	
	/***************USAGE***************
     oscgen outfile dur srate num_chans, amp, freq, wavetype, num_oscs
     **********************************/
    if(argc < OSC_GEN_NARGS)
    {
        fprintf(stderr, "insufficient arguments.\n"
				"usage: oscgen outfile dur srate num_chans amp freq wavetype   num_oscs\n"
				"where wavetype =:\n"
				"			0 = square\n"
				"			1 = triangle\n"
				"			2 = sawtooth up\n"
				"			3 = sawtooth down\n"
				"			4 = pulse\n"
				"dur		= duration of outfile (seconds)\n"
				"srate		= required sample rate of outfile\n"
				"amp		= amplitude value or breakpoint file (0 < amp <= 1.0)\n"
				"freq		= frequency value (freq > 0) or breakpoint file.\n"
				"num_oscs	= number of oscillators (where 1 is normal sine wave)\n"
				);
        exit(1);
    }
	
	/********* initialize all dynamic resources to default states **********/
	PSF_PROPS outprops;
	int ofd = -1;
	int error = 0;
	PSF_CHPEAK* peaks = NULL;	
	psf_format outformat =  PSF_FMT_UNKNOWN;
	long nframes = NFRAMES;
	float* outframe = NULL;
	double amp,freq,dur;
	long outframes,nbufs = 0,num_oscs,i; //num_oscs holds number of oscillators for additive synthesis
	long remainder;
	int samptype = PSF_SAMP_16;
	OSCIL** oscs = NULL; //an array of oscil pointers (because additive synthesis requires more than one oscillators)
	double* oscamps = NULL; //for oscillator bank amplitudes and frequencies
	double* oscfreqs = NULL;
	double phase = 0.0;
	double ampfac,freqfac,ampadjust; //Because we cannot know the number of oscillators before hand, we need to calculate it cumulatively
	//Recall that adding sine waves amplitude is as simple as numerically adding the amplitude values. So 1 + 3 + 5 + 7 will yield an amp
	//value of 16. However, we want amp to be between -1 and 1 --> therefore total amp is 1/(1^2) + 1/(3^2) + 1/(5^2) + 1/(7^2) = 1 + 1/9 + 1/25 + 1/49 = 1.172
	//General Rule : totalamp = Sum (from n = 0 to N) 1/(2n-1)^2 where N = number of oscillators
	
	FILE* fpfreq = NULL;
	BRKSTREAM *freqstream = NULL, *ampstream = NULL;
	unsigned long brkfreqSize = 0;
	double minval = 0.0,maxval= 0.0;
	FILE* fpamp = NULL;
	unsigned long brkampSize = 0;
	double amp_minval = 0.0,amp_maxval= 0.0;
	clock_t starttime, endtime;
	int wavetype;
	oscs = NULL;
	
	
	/********** error checking on wave type **********/
    wavetype = atoi(argv[OSC_GEN_WAVETYPE]);
    if(wavetype < OSC_GEN_SQUARE || wavetype > OSC_GEN_PULSE)
    {
        fprintf(stderr, "Error: wavetype not defined\n");
        goto exit;
    }
	
	/********** error checking on number of oscillators **********/
	num_oscs = atoi(argv[OSC_GEN_NUM_OSCS]);
	if(num_oscs <= 0)
	{
		fprintf(stderr, "Error: number of oscillators must be positive\n");
		goto exit;
	}
	
	/********** define outfile format (with embedded error checking) **********/
	outprops.srate = atoi(argv[OSC_GEN_SRATE]);
	if(outprops.srate <= 0)
	{
		fprintf(stderr, "Error: srate must be positive\n");
		goto exit;
	}
	outprops.chans = atoi(argv[OSC_GEN_NUM_CHANS]);
	if(outprops.chans <= 0)
	{
		fprintf(stderr, "Error: number of channels must be positive\n");
		goto exit;
	}
	outprops.samptype = (psf_stype) samptype;
	outprops.chformat = STDWAVE;
	
	dur = atof(argv[OSC_GEN_DUR]);
    if(dur <= 0.0)
    {
        fprintf(stderr, "Error: duration must be positive (>=0)\n");
        goto exit;
    }
    outframes = (unsigned long) (dur * outprops.srate + 0.5);
    nbufs = outframes / nframes;
    remainder = outframes - nbufs * nframes;
    if(remainder > 0) nbufs++;
	
	
	/********** open breakpoint files, or set constants if no breakpoint file is defined **********/
    /***** amp (amplitude) breakpoints file *****/
    fpamp = fopen(argv[OSC_GEN_AMP], "r"); //try opening specified argument
    if(fpamp == NULL)
    {
        //meaning that either argument OSC_GEN_AMP is not a breakpoint file or error opening
        amp = atof(argv[OSC_GEN_AMP]);
        if(amp <= 0.0 || amp > 1)
        {
            fprintf(stderr, "Error: amplitude must be between 0 and 1 (0 <= amp <= 1)\n");
            goto exit;
        }
    }
    else
    {
        //meaning that argument OSC_GEN_AMP is indeed a breakpoint file
        ampstream = bps_newstream(fpamp,outprops.srate,&brkampSize);
        if(ampstream == NULL)
        {
            fprintf(stderr, "Error reading amplitude breakpoint file: %s\n", argv[OSC_GEN_AMP]);
            goto exit;
        }
        if(bps_getminmax(ampstream, &minval, &maxval))
        {
            //meaning we called upon said function that will put into minval and maxval appropriate values
            //if said function returns 1 --> we have error
            fprintf(stderr, "Error reading range breakpoint file: %s\n", argv[OSC_GEN_AMP]);
            goto exit;
        }
        if(minval < 0.0 || minval > 1.0 || maxval < 0.0 || maxval > 1.0)
        {
            fprintf(stderr, "Error: amplitude values out of range in breakpoint file\n");
            goto exit;
        }
    }
    
    /***** freq (frequency) breakpoints file *****/
    fpfreq = fopen(argv[OSC_GEN_FREQ], "r");
    if(fpfreq == NULL)
    {
        //meaning that either argument OSC_GEN_FREQ is not a breakpoint file or error opening
        freq = atof(argv[OSC_GEN_FREQ]);
        if(freq <= 0.0)
        {
            fprintf(stderr, "Error: frequency must be positive\n");
            goto exit;
        }
    }
    else
    {
        //meaning that argument OSC_GEN_FREQ is indeed a breakpoint file
        freqstream = bps_newstream(fpfreq, outprops.srate, &brkfreqSize);
        if(freqstream == NULL)
        {
            fprintf(stderr, "Error reading frequncy from breakpoint file: %s\n", argv[OSC_GEN_FREQ]);
            goto exit;
        }
        if(bps_getminmax(freqstream, &minval, &maxval))
        {
            fprintf(stderr, "Error reading range from breakpoint file: %s\n", argv[OSC_GEN_FREQ]);
            goto exit;
        }
        if(minval <= 0.0)
        {
            fprintf(stderr, "Error in breakpoing file; Frequency values must be positive\n");
            goto exit;
        }
		
    }
	
	/********** Initialize portsf **********/
    if(psf_init())
    {
        fprintf(stderr, "Error: Unable to initialize portsf\n");
        goto exit;
    }
    
    /********** Sample Frame Buffer Memory Aloocation **********/
    outframe = (float*)malloc(nframes * outprops.chans * sizeof(float));
    if(outframe == NULL)
    {
        fprintf(stderr, "Error: No memory for outframe buffer\n");
        goto exit;
    }
    
    /********** Check for correct output extension **********/
    outformat = psf_getFormatExt(argv[OSC_GEN_OUTFILE]);
    if(outformat == PSF_FMT_UNKNOWN)
    {
        fprintf(stderr, "Error: invalid extension\n"
                "Possible extensions:\t *.wav\t *.aiff\t *.aif\t *.afc\t *.aifc\n");
        goto exit;
    }
    //reach here meaning outformat is viable extension
    outprops.format = outformat;
	
	/********** create oscillator amp and oscillator freq arrays **********/
	oscamps = (double*)malloc(num_oscs * sizeof(double));
	oscfreqs = (double*)malloc(num_oscs * sizeof(double));
	if(!oscamps || !oscfreqs)
	{
		fprintf(stderr, "Error: Failed memory allocation for oscillator amp/freq arrays\n");
		goto exit;
	}
	
	/********** create amp/freq data for requested waveshape, normalized to 1 **********/
	ampfac = 1.0;
	freqfac = 1.0;
	ampadjust = 0.0;
	
	switch(wavetype)
	{
		case(OSC_GEN_SQUARE):
			for(i=0; i<num_oscs; i++)
			{
				oscamps[i] = ampfac;
				oscfreqs[i] = freqfac;
				freqfac += 2.0;
				ampadjust += ampfac;
				ampfac = 1.0 / freqfac;
			}
			//normalize
			for(i=0; i<num_oscs; i++)
				oscamps[i] /= ampadjust;
			break;
		case(OSC_GEN_TRIANGLE):
			for(i=0; i<num_oscs; i++)
			{
				oscamps[i] = ampfac;
				oscfreqs[i] = freqfac;
				freqfac += 2.0;
				ampadjust += ampfac;
				ampfac = 1.0 / (freqfac * freqfac);
			}
			//normalize
			for(i=0; i<num_oscs; i++)
				oscamps[i] /= ampadjust;
			phase = 0.25;
			break;
		case(OSC_GEN_SAW_UP):
		case(OSC_GEN_SAW_DOWN):
			for(i=0; i<num_oscs; i++)
			{
				oscamps[i] = ampfac;
				oscfreqs[i] = freqfac;
				freqfac += 1.0;
				ampadjust += ampfac;
				ampfac = 1.0 / freqfac;
			}
			if(wavetype == OSC_GEN_SAW_UP)
				ampadjust = -ampadjust;
			//normalize
			for(i=0; i<num_oscs; i++)
				oscamps[i] /= ampadjust;
			break;
		case(OSC_GEN_PULSE):
			freqfac = 1.0;
			ampfac = 1.0 / num_oscs;
			ampadjust = 0.0;
			for(i=0; i<num_oscs; i++)
			{
				oscamps[i] = ampfac;
				oscfreqs[i] = freqfac;
				freqfac += 1.0;
			}
			break;			
	}
	
	/********** Create Oscillator Bank **********/
	oscs = (OSCIL**)malloc(num_oscs * sizeof(OSCIL*));
	if(!oscs)
	{
		fprintf(stderr, "Error: Failed memory allocation on creation of oscillator bank\n");
		goto exit;
	}
	for(i=0; i< num_oscs; i++)
	{
		oscs[i] = new_oscil_with_phase(outprops.srate, phase);
		if(!oscs[i])
		{
			fprintf(stderr, "Error: Failed memory allocation on creation of oscillator object\n");
			goto exit;
		}
	}
	
	
	/********** Initialize Peak nfo **********/
    peaks = (PSF_CHPEAK*)malloc(outprops.chans * sizeof(PSF_CHPEAK));
    if(peaks == NULL)
    {
        fprintf(stderr, "Error: No memory for peaks\n");
        goto exit;
    }
    
    /********** Create and Handle outfile **********/
    ofd = psf_sndCreate(argv[OSC_GEN_OUTFILE],&outprops, 0, 0, PSF_CREATE_RDWR); //refer to sf2float for documentation
    if(ofd < 0)
    {
        fprintf(stderr, "Error: Unable to create outfile\n");
        goto exit;
    }
    
    /********** Audio Processing **********/
	printf("Reached here\t\t numbufs: %lu; numframes: %lu\n", nbufs, nframes);
	starttime = clock();
	for(i=0; i<nbufs; i++)
	{
		long j,k,l;
		double val;
		if(i == nbufs -1)
			nframes = remainder;
		
		for(j=0; j< nframes; j++)
		{
			if(freqstream)
				freq = bps_tick(freqstream);
			if(ampstream)
				amp = bps_tick(freqstream);
			val = 0.0;
			for(l=0; l<num_oscs; l++)
			{
				val += oscamps[l] * sine_wave_node(oscs[l], freq * oscfreqs[l]);
			}
			for(k=0; k<outprops.chans; k++)
			{
				outframe[j*outprops.chans +k] = (float)(val * amp);
			}
		} //end j for-loop
		if(psf_sndWriteFloatFrames(ofd, outframe, nframes) != nframes)
		{
			fprintf(stderr, "Error writing to outfile");
			break;
		}
	}//end i for-loop
	endtime = clock();
	
exit:
	fprintf(stdout, "Task Finished \nElapsed Time: %f seconds\n", (endtime - starttime)/(double)CLOCKS_PER_SEC);
    /***** Report PEAK values to user *****/
    if(psf_sndReadPeaks(ofd,peaks,NULL) > 0){
		long i;
		double peaktime;
		printf("PEAK information:\n");	
        for(i=0;i < outprops.chans;i++){
            peaktime = (double) peaks[i].pos / (double) outprops.srate;
            printf("CH %d:\t%.4f at %.4f secs\n", (int)(i+1), peaks[i].val, peaktime);
        }
	}

	if(ofd >= 0)
		psf_sndClose(ofd);
	if(outframe)
		free(outframe);
	if(peaks)
		free(peaks);
	
	if(oscs){
		for(i=0;i < num_oscs;i++){
			free(oscs[i]);
		}
		free(oscs);
	}
	if(oscamps)
		free(oscamps);
	if(oscfreqs)
		free(oscfreqs);
	psf_finish();
	exit(1);
}

void summation_function(int argc, char *argv[])
{
	//This function is to compute the sum of the series a^n for some a>1 and n
	//It is for CS473 SP12 Quiz 3

	int a = 2;
	int n = 25;
	
	int i = 0; //index
	
	int a_exp_n = 0;
	
	for(i=0; i<=n; i++)
	{
		int additive = pow(a, i);
		printf("Additive: %d^%d = %d\n", a, i, additive);
		a_exp_n += additive;
		printf("Total at %d: %d\n", i, a_exp_n);
	}
	
	//Here we test to see which equation approximates the sum
	int option_A = (pow(a, n+1) + 1)/(a-1); // (a^(n+1) +1) / (a-1)
	int option_B = (pow(a, n+1) - 1)/(a-1); // (a^(n+1) -1) / (a-1)
	int option_C = (pow(a, n+1) - 1)/(a+1); // (a^(n+1) -1) / (a+1)
	int option_D = pow(a, n+1); // a^(n+1)
	
	printf("\n");
	printf("Option A: %d\n", option_A);
	printf("Option B: %d\n", option_B);
	printf("Option C: %d\n", option_C);
	printf("Option D: %d\n", option_D);
	
	printf("\n");
	if(a_exp_n == option_A) printf("Option_A is the correct approximation\n");
	else if(a_exp_n == option_B) printf("Option_B is the correct approximation\n");
	else if(a_exp_n == option_C) printf("Option_C is the correct approximation\n");
	else if(a_exp_n == option_D) printf("Option_D is the correct approximation\n");
	
}

void secant_method(int argc, char *argv[])
{
	//Secant method is x_n = x_n-1 - f(x_n-1)* ((x_n-1 - x_n-2) / (f(x_n-1) - f(x_n-2)))
	//From quiz 3, x1 = 2.1; x0 = 2.2 where x_n-1 = x1 and x_n-2 = x0
	//From quiz 3, another problem, x0 = 1.0 and x1 = 1.1
	
	float x1 = 1.0;
	float x0 = 1.1;
	float result = x1 - fnx(x1) * ( (x1 - x0) / (fnx(x1) - fnx(x0)));
	printf("The result is: %f\n", result);
}

float fnx(float x)
{
	//the original function was f(x) = x^2 + x -6 for one Q
	//the original function was f(x) = x^2 - 1.1 for another Q
	float y = pow(x,2) - 1.1;
	return y;
}


void oddSumTest(int argc, char *argv[])
{
	printf("%s\n", argv[1]);
	int num = atoi(argv[1]);
	int sum = 0;
	
	if (num%2 == 0) num = num -1;
	
	int i;
	for(i = 1; i <= num - 2; i=i+2)
	{
		sum += i;
		printf("%d + ", i);
	}
	printf("%d =  %d\n", num, sum + num);
	
}

void pointerTest()
{
	int *p, q, x;
	x=10;
	p=&x;
	*p = x+1;
	q=x;
	printf("q = %d\n", q);
	
	int ptr[2];
	ptr[1] = 1;
	ptr[2] = 2;
	printf("%d\n", ptr[1]);
	printf("%d\n", ptr[2]);
	
	char str[5];
	strcpy (str, "Hello");
	printf("str: %s\n", str);
	
	int x1 = 10;
	int y1, y2;
	y1 = x1++; 
	printf("%d\n", y1);
	
	int x2 = 10;
	y2 = ++x2;
	printf("%d\n", y2);
	
	int *p2 = (int *)500;
	printf("%lu\n", (size_t)p2);

}

void aQRec()
{
	int bufferSize;
	for(bufferSize=0; bufferSize < 5; ++bufferSize)
		printf("++bufferSize value: %d\n", bufferSize);
	printf("++After: %d\n", bufferSize);
	for (bufferSize=0; bufferSize < 5; bufferSize++) 
		printf("bufferSize++ value: %d\n", bufferSize);
	printf("After++: %d\n", bufferSize);
	
	printf("Currently Working. Press <return> to exit.\n");
	getchar(); //Note that on default, it will close on <return>
	
}

float m_POW(int a, int b)
{
	if(b==0) return 1;
	else if(b==1) return a;
	int i, temp = a;
	for(i=0; i<=b-2; i++)
		temp = a * temp;
	return temp;
}













