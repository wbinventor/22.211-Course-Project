
/*-------------------------------------------------------------------------*/
/**
  @file		gnuplot_i.h
  @author	N. Devillard
  @date		Sep 1998
  @version	$Revision: 1.11 $
  @brief	C interface to gnuplot.

  gnuplot is a freely available, command-driven graphical display tool for
  Unix. It compiles and works quite well on a number of Unix flavours as
  well as other operating systems. The following module enables sending
  display requests to gnuplot through simple C calls.
  
*/
/*--------------------------------------------------------------------------*/

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <unistd.h>
#include <stdarg.h>
#include <string>


/** Maximal number of simultaneous temporary files */
#define GP_MAX_TMP_FILES    64
/** Maximal size of a temporary file name */
#define GP_TMP_NAME_SIZE    512


/**
 *@typedef	gnuplot_ctrl
 *@brief	gnuplot session handle (opaque type).

  This structure holds all necessary information to talk to a gnuplot
  session. It is built and returned by gnuplot_init() and later used
  by all functions in this module to communicate with the session, then
  meant to be closed by gnuplot_close().

  This structure is meant to remain opaque, you normally do not need
  to know what is contained in there.
 */
typedef struct _GNUPLOT_CTRL_ {

	/** Pipe to gnuplot process */
    FILE    * gnucmd ;
    
    /** Number of currently active plots */
    int       nplots ;

    /** Current plotting style */
    char      pstyle[32] ;

    /** Name of temporary files */
    char      to_delete[GP_MAX_TMP_FILES][GP_TMP_NAME_SIZE] ;

    /** Number of temporary files */
    int       ntmp ;
} gnuplot_ctrl ;

char * gnuplot_get_program_path(char* pname);
gnuplot_ctrl * gnuplot_init(void);
void gnuplot_close(gnuplot_ctrl* handle);
void gnuplot_cmd(gnuplot_ctrl*  handle, char*  cmd, ...);
void gnuplot_setstyle(gnuplot_ctrl* h, char* plot_style);
void gnuplot_set_xlabel(gnuplot_ctrl* h, char* label);
void gnuplot_set_ylabel(gnuplot_ctrl* h, char* label);
void gnuplot_set_yrange(gnuplot_ctrl* h, int lower, int upper);
void gnuplot_set_yrange(gnuplot_ctrl* h, float lower, float upper);
void gnuplot_set_yrange(gnuplot_ctrl* h, double lower, double upper);
void gnuplot_set_xrange(gnuplot_ctrl* h, int lower, int upper);
void gnuplot_set_xrange(gnuplot_ctrl* h, float lower, float upper);
void gnuplot_set_xrange(gnuplot_ctrl* h, double lower, double upper);
void gnuplot_resetplot(gnuplot_ctrl* h);
void gnuplot_saveplot(gnuplot_ctrl* h, char* filename);

void gnuplot_plot_x(gnuplot_ctrl* handle, int* d, int n, char* title);
void gnuplot_plot_x(gnuplot_ctrl* handle, float* d, int n, char* title);
void gnuplot_plot_x(gnuplot_ctrl* handle, double* d, int n, char* title);

void gnuplot_plot_xy(gnuplot_ctrl* handle, int* x, int* y,
											int n, char* title);
void gnuplot_plot_xy(gnuplot_ctrl* handle, float* x, float* y,
											int n, char* title);
void gnuplot_plot_xy(gnuplot_ctrl* handle, double* x, double* y,
											int n, char* title);

void gnuplot_plot_once(char* title, char* style, char* label_x, char* label_y,
													float* x, float* y, int n);
void gnuplot_plot_once(char* title, char* style, char* label_x, char* label_y,
													float* x, float* y, int n);
void gnuplot_plot_once(char* title, char* style, char* label_x, char* label_y,
													float* x, float* y, int n);

void gnuplot_plot_slope(gnuplot_ctrl* handle, float a, float b, char* title);

void gnuplot_plot_equation(gnuplot_ctrl* h, char* equation, char* title) ;

void gnuplot_plot_xyz(gnuplot_ctrl* handle, int* x, int* y, int* z,
														int n, char* title);
void gnuplot_plot_xyz(gnuplot_ctrl* handle, float* x, float* y, float* z,
														int n, char* title);
void gnuplot_plot_xyz(gnuplot_ctrl* handle, double* x, double* y, double* z,
														int n, char* title);

#endif
