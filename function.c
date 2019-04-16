#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ----------------------------------------------------------------------
** Test input signal contains 1000Hz + 15000 Hz
** ------------------------------------------------------------------- */

double signalMean(double * sigSourceArray, int sigLength)
{
    double mean = 0.0;
    for(int i = 0; i < sigLength; i++)
    {
        mean = mean + sigSourceArray[i];
    }
    mean = mean / (double)sigLength;

    return mean;
}

double signalVariance(double * sigSourceArray, int sigLength, double sigMean)
{
    double variance = 0.0;

     for(int i = 0; i < sigLength; i++)
    {
        variance += pow((sigSourceArray[i] - sigMean), 2);
    }

    variance /= (sigLength -1);

    return variance;
}

double signalStandardDeviation(double sigVariance)
{
    double std = sqrt(sigVariance);
    return std;
}

void signalConvolution(double* sigSourceArray, int sigLength, double* impResponseArray, int impLength, double* sigOutputArray)
{
    for(int i = 0; i < (sigLength + impLength); i++)
    {
        sigOutputArray[i] = 0;
    }

    for(int i = 0; i < sigLength; i++)
    {
        for (int j = 0; j < impLength; j++)
        {
            sigOutputArray[i + j] += sigSourceArray[i] * impResponseArray[j];
        }
    }
}

void runningSum(double* sigSourceArray, int sigLength, double* sigOutputArray)
{
    for(int i = 0; i < sigLength; i++)
    {
        sigOutputArray[i] = sigOutputArray[i - 1] + sigSourceArray[i];
    }
}

void dft(double* sigSourceArray, int sigLength, double* sigOutputArrayReal, double* sigOutputArrayImm)
{
    int i,j,k;
    double PI = 3.14159265359;

    for(j = 0; j < sigLength / 2; j++)
    {
        sigOutputArrayReal[j] = 0;
        sigOutputArrayImm[j] = 0;
    }

    for(k = 0; k < sigLength/2; k++)
    {
        for(i = 0; i < sigLength; i++)
        {
            sigOutputArrayReal[k] += sigSourceArray[i] * cos(2 * PI * k * i / sigLength);
            sigOutputArrayImm[k] -= sigSourceArray[i] * sin(2 * PI * k * i / sigLength);
        }
    }
}

void inverse_dft(double* sigOutputArray, int sigLength, double* sigSourceArrayReal, double* sigSourceArrayImm)
{
    double PI = 3.14159265359;

    int i, k;

    for(k = 0; k < sigLength / 2; k++)
    {
        sigSourceArrayReal[k] = sigSourceArrayReal[k] / (sigLength / 2);
        sigSourceArrayImm[k] = -sigSourceArrayImm[k] / (sigLength / 2);
    }

    sigSourceArrayReal[0] = sigSourceArrayReal[0] / 2;
    sigSourceArrayImm[0] = -sigSourceArrayImm[0] /  2;

    for(i = 0; i < sigLength; i++)
    {
        sigOutputArray[i] = 0;
    }

    for(k = 0; k < (sigLength / 2); k++)
    {
        for(i = 0; i < sigLength; i++)
        {
            sigOutputArray[i] += (sigSourceArrayReal[i] * cos((2 * PI * k * i) / sigLength));
            sigOutputArray[i] += (sigSourceArrayImm[i] * sin((2 * PI * k * i) / sigLength));
        }
    }
}

void complex_dft(double* sigSourceArrayReal_timeDomain, double* sigSourceArrayImm_timeDomain, double* sigOutputArrayReal_freqDomain, double* sigOutputArrayImm_freqDomain, int sigLength)
{
    double PI = 3.14159265359;
    double SR, SI;

    for(int k = 0; k < sigLength; k++)
    {
        for(int i = 0; i < sigLength - 1; i++)
        {
            SR = cos((2 * PI * k * i) / sigLength);
            SI = -sin((2 * PI * k * i) / sigLength);

            sigOutputArrayReal_freqDomain[k] += (sigSourceArrayReal_timeDomain[i] * SR - sigSourceArrayImm_timeDomain[i] * SI);
            sigOutputArrayImm_freqDomain[k] += (sigSourceArrayImm_timeDomain[i] * SI - sigSourceArrayImm_timeDomain[i] * SR);
        }
    }
}

/**
    note: cutoff freq is normalized, must be between 0 and 0.5 (represent Nyquist freq.)
    the input signal in this example was sampled at 48 kHz, so Nyquist freq. is 24 kHz,
    to understand:
    - if
    24 kHz ----> 0.5
    - so
    10 kHz ----> 0.2
    because (10/24) x 0.5 = 0.2
**/
void lowpass_windowed_sinc_ftr(double* sig_sourceArray, double* sig_outputArray, double* fltr_kernel_outputArray,
                                double cutoff_freq, int filter_length, int input_signal_length)
{
    double PI = 3.14159265359;

    // calculate the low pass filter kernel:
    for(int i = 0; i < filter_length; i++)
    {
        if(i - filter_length / 2 == 0)
        {
            fltr_kernel_outputArray[i] = 2 * PI * cutoff_freq;
        }
        if(i - filter_length / 2 != 0)
        {
            fltr_kernel_outputArray[i] = sin(2 * PI * cutoff_freq * (i - (filter_length / 2))) / (i - (filter_length / 2));
            fltr_kernel_outputArray[i] = fltr_kernel_outputArray[i] * (0.54 - 0.46 * cos(2 * PI * i / filter_length)); // 0.54 and 0.46 are specific parameter of Hamming filter
        }
    }

    // convolve input signal with filter kernel
    for(int j = filter_length; j < input_signal_length;  j++)
    {
        sig_outputArray[j] = 0;
        for(int i = 0; i < filter_length; i ++)
        {
            sig_outputArray[j] = sig_outputArray[j] + (sig_sourceArray[j - i] * fltr_kernel_outputArray[i]);
        }
    }
}

/**
    0,1 kHz ----> 0.002
    because (0,1/24) x 0.5 = 0.002
    5 kHz ----> 0.11
    because (5/24) x 0.5 = 0.11
**/
void bandpass_windowed_sinc_ftr (double* sig_sourceArray, double* sig_outputArray, double* fltr_kernel_outputArray,
                                 double low_cutoff_freq, double high_cutoff_freq, int filter_length,
                                 int input_signal_length, double* low_cutoff_buffer, double*  high_cutoff_buffer)
{
    double PI = 3.14159265359;

    // calculate the low pass filter kernel:
    for(int i = 0; i < filter_length; i++)
    {
        if(i - filter_length / 2 == 0)
        {
            low_cutoff_buffer[i] = 2 * PI * low_cutoff_freq;
        }

        if(i - filter_length / 2 != 0)
        {
            low_cutoff_buffer[i] = sin(2 * PI * low_cutoff_freq * (i - filter_length / 2)) / (i - filter_length / 2);
            low_cutoff_buffer[i] = low_cutoff_buffer[i] * ((0.42 - 0.5 * cos(2 * PI * i / filter_length)) + 0.08 * cos(4 * PI * i / filter_length));
        }
    }

    // calculate the high pass filter kernel, by a spectral inverted low pass:
    for(int i = 0; i < filter_length; i++)
    {
        if(i - filter_length / 2 == 0)
        {
            high_cutoff_buffer[i] = 2 * PI * high_cutoff_freq;
        }

        if(i - filter_length / 2 != 0)
        {
            high_cutoff_buffer[i] = sin(2 * PI * high_cutoff_freq * (i - filter_length/2)) / (i - filter_length / 2);
            high_cutoff_buffer[i] = high_cutoff_buffer[i] * (0.42 - 0.5 * cos(2 * PI * i / filter_length) + 0.08 * cos(4 * PI * i / filter_length));
        }
    }

    for(int i = 0; i < filter_length; i++)
    {
        high_cutoff_buffer[i] = -high_cutoff_buffer[i];
    }

    high_cutoff_buffer[filter_length / 2] = high_cutoff_buffer[filter_length / 2] + 1;

    // create kernel for band pass filter by a spectral inverted band reject filter:
    for(int i = 0; i < filter_length; i++)
    {
        fltr_kernel_outputArray[i] = low_cutoff_buffer[i] +  high_cutoff_buffer[i];
    }

    for(int i = 0; i < filter_length; i++)
    {
         fltr_kernel_outputArray[i] = -fltr_kernel_outputArray[i];
    }
    fltr_kernel_outputArray[filter_length / 2] = fltr_kernel_outputArray[filter_length / 2] + 1;

    // now convolve input signal with band pass filter kernel
    for(int j = filter_length; j < input_signal_length;  j++)
    {
        sig_outputArray[j] = 0;
        for(int i = 0; i < filter_length; i ++)
        {
            sig_outputArray[j] = sig_outputArray[j] + (sig_sourceArray[j - i] * fltr_kernel_outputArray[i]);
        }
    }
}
