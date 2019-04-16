#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "waveform.h"
#include "function.c"

#define SIG_LENGTH          320
#define SIG_2_LENGTH        501

#define IMP_RSP_LENGTH      29

#define KERNEL_LENGTH       29

// define variables for statistic parameters
double _mean;
double variance;
double standardDeviation;

// define variables for convolution
double outputSignal[SIG_LENGTH + IMP_RSP_LENGTH];
double outputSignal2[SIG_LENGTH];

// define variables for Fourier Analysis
double outputSignalReal[SIG_LENGTH/2];
double outputSignalImm[SIG_LENGTH/2];

double output_idft[SIG_LENGTH];

double output_complex_dft_Real[SIG_LENGTH];
double output_complex_dft_Imm[SIG_LENGTH];

// define variables for digital filter design
double outputSignal_fromLowFilter[SIG_LENGTH - KERNEL_LENGTH];
double outputFilterKernel[KERNEL_LENGTH];

double outputSignal_fromBandpassFilter[SIG_LENGTH - KERNEL_LENGTH];
double outputBandpassFilterKernel[KERNEL_LENGTH];
double low_cutoff_buffer[KERNEL_LENGTH];
double high_cutoff_buffer[KERNEL_LENGTH];



int main()
{
    // Mean
    _mean = signalMean(&InputSignal_f32_1kHz_15kHz[0], SIG_LENGTH);

    printf(" \nmean = \n%f", _mean);

    // Variance
    variance = signalVariance(&InputSignal_f32_1kHz_15kHz[0], SIG_LENGTH, _mean);

    printf(" \nvariance = \n%f", variance);

    // Standard Deviation
    standardDeviation =  signalStandardDeviation(variance);

    printf(" \nstd = \n%f", standardDeviation);


    FILE* input_sig_fptr; // create files to plot this signal with gnuplot
    FILE* imp_rsp_fptr;

    input_sig_fptr = fopen("input_signal.dat", "w");
    imp_rsp_fptr = fopen("impulse_response.dat", "w");

    for(int i = 0; i < SIG_LENGTH; i++)
    {
        fprintf(input_sig_fptr, "\n%f", InputSignal_f32_1kHz_15kHz[i]);
    }

    fclose(input_sig_fptr);

    for(int i = 0; i < IMP_RSP_LENGTH; i++)
    {
        fprintf(imp_rsp_fptr, "\n%f", Impulse_response[i]);
    }

    fclose(imp_rsp_fptr);

    // Convolution
    signalConvolution (&InputSignal_f32_1kHz_15kHz[0], SIG_LENGTH, &Impulse_response[0], IMP_RSP_LENGTH, &outputSignal[0]);

    FILE* output_sig_fptr;

    output_sig_fptr = fopen("output_signal.dat", "w");

    for(int i = 0; i < SIG_LENGTH + IMP_RSP_LENGTH; i++)
    {
        fprintf(output_sig_fptr, "\n%f", outputSignal[i]);
    }

    fclose(output_sig_fptr);


    FILE* output_sig_fptr2;


    runningSum(&InputSignal_f32_1kHz_15kHz[0], SIG_LENGTH, &outputSignal2[0]);

    output_sig_fptr2 = fopen("output_signal2.dat", "w");

    for(int i = 0; i < SIG_LENGTH; i++)
    {
        fprintf(output_sig_fptr2, "\n%f", outputSignal2[i]);
    }

    fclose(output_sig_fptr2);

    // Discrete Fourier Transform
    FILE* fptr;
    FILE* fptr2;
    FILE* fptr3;

    dft(&InputSignal_f32_1kHz_15kHz[0], SIG_LENGTH, &outputSignalReal[0], &outputSignalImm[0]);

    fptr = fopen("input_signal.dat", "w");
    fptr2 = fopen("output_signalReal.dat", "w");
    fptr3 = fopen("output_signalImm.dat", "w");

    for(int i = 0; i < SIG_LENGTH; i++)
    {
        fprintf(fptr, "\n%f", InputSignal_f32_1kHz_15kHz[i]);
    }

    for(int i = 0; i < SIG_LENGTH / 2; i++)
    {
        fprintf(fptr2, "\n%f", outputSignalReal[i]);
        fprintf(fptr3, "\n%f", outputSignalImm[i]);
    }

    fclose(fptr);
    fclose(fptr2);
    fclose(fptr3);

    // Inverser dft
    FILE* fptr4;

    inverse_dft(&output_idft[0], SIG_LENGTH, &outputSignalReal[0], &outputSignalImm[0]);

    fptr4 = fopen("output_inversedft.dat", "w");

    for(int i = 0; i < SIG_LENGTH; i++)
    {
        fprintf(fptr4, "\n%f", output_idft[i]);
    }

    fclose(fptr4);

    // Complex dft
    FILE* fptr5;
    FILE* fptr6;
    FILE* fptr7;
    FILE* fptr8;

    complex_dft(&InputSignal_501pts_20Hz_sigReal[0], &InputSignal_501pts_20Hz_sigImm[0], &output_complex_dft_Real[0], &output_complex_dft_Imm[0], SIG_2_LENGTH);

    fptr5 = fopen("input_signal_complexdft_Real.dat", "w");
    fptr6 = fopen("input_signal_complexdft_Imm.dat", "w");
    fptr7 = fopen("output_complexdft_Real.dat", "w");
    fptr8 = fopen("output_complexdft_Imm.dat", "w");

    for(int i = 0; i < SIG_2_LENGTH; i++)
    {
        fprintf(fptr5, "\n%f", InputSignal_501pts_20Hz_sigReal[i]);
        fprintf(fptr6, "\n%f", InputSignal_501pts_20Hz_sigImm[i]);
        fprintf(fptr7, "\n%f", output_complex_dft_Real[i]);
        fprintf(fptr8, "\n%f", output_complex_dft_Imm[i]);
    }

    fclose(fptr5);
    fclose(fptr6);
    fclose(fptr7);
    fclose(fptr8);

    // Digital filter
    FILE* fptr9;
    FILE* fptr10;
    FILE* fptr11;

    lowpass_windowed_sinc_ftr(&InputSignal_f32_1kHz_15kHz[0], &outputSignal_fromLowFilter[0], &outputFilterKernel[0], 0.2, KERNEL_LENGTH, SIG_LENGTH);

    fptr9 = fopen("output_lowFilter.dat", "w");
    fptr10 = fopen("input_lowFilter.dat", "w");
    fptr11 = fopen("output_kernel.dat", "w");

    for(int i = 0; i < SIG_LENGTH; i++)
    {
        fprintf(fptr10, "\n%f", InputSignal_f32_1kHz_15kHz[i]);
        if(i >= KERNEL_LENGTH)
        {
            fprintf(fptr9, "\n%f", outputSignal_fromLowFilter[i]);
        }
        if(i <= KERNEL_LENGTH)
        {
            fprintf(fptr11, "\n%f", outputFilterKernel[i]);
        }

    }

    fclose(fptr9);
    fclose(fptr10);
    fclose(fptr11);

    FILE* fptr12;
    FILE* fptr13;
    FILE* fptr14;

    bandpass_windowed_sinc_ftr(&InputSignal_f32_1kHz_15kHz[0], &outputSignal_fromBandpassFilter[0], &outputBandpassFilterKernel, 0.002, 0.11, KERNEL_LENGTH, SIG_LENGTH, &low_cutoff_buffer, &high_cutoff_buffer);

    fptr12 = fopen("output_bandpassFilter.dat", "w");
    fptr13 = fopen("input_bandpassFilter.dat", "w");
    fptr14 = fopen("output_bandpasskernel.dat", "w");

    for(int i = 0; i < SIG_LENGTH; i++)
    {
        fprintf(fptr13, "\n%f", InputSignal_f32_1kHz_15kHz[i]);
        if(i >= KERNEL_LENGTH)
        {
            fprintf(fptr12, "\n%f", outputSignal_fromBandpassFilter[i]);
        }
        if(i <= KERNEL_LENGTH)
        {
            fprintf(fptr14, "\n%f", outputBandpassFilterKernel[i]);
        }
    }

    fclose(fptr12);
    fclose(fptr13);
    fclose(fptr14);


    return 0;
}
