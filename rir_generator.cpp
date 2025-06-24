#define _USE_MATH_DEFINES

#include "matrix.h"
#include "mex.h"
#include "math.h"
#include "rir_generator_core.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs == 0)
    {
        mexPrintf("--------------------------------------------------------------------\n"
            "| Room Impulse Response Generator                                  |\n"
            "|                                                                  |\n"
            "| Computes the response of an acoustic source to one or more       |\n"
            "| microphones in a reverberant room using the image method [1,2].  |\n"
            "|                                                                  |\n"
            "| Author    : dr.ir. Emanuel Habets (e.habets@ieee.org)            |\n"
            "|                                                                  |\n"
            "| Version   : 2.3.20241022                                         |\n"
            "|                                                                  |\n"
            "| Copyright (C) 2003-2025 E.A.P. Habets                            |\n"
            "|                                                                  |\n"
            "| [1] J.B. Allen and D.A. Berkley,                                 |\n"
            "|     Image method for efficiently simulating small-room acoustics,|\n"
            "|     Journal Acoustic Society of America,                         |\n"
            "|     65(4), April 1979, p 943.                                    |\n"
            "|                                                                  |\n"
            "| [2] P.M. Peterson,                                               |\n"
            "|     Simulating the response of multiple microphones to a single  |\n"
            "|     acoustic source in a reverberant room, Journal Acoustic      |\n"
            "|     Society of America, 80(5), November 1986.                    |\n"
            "--------------------------------------------------------------------\n\n"
            "function [h, beta_hat] = rir_generator(c, fs, r, s, L, beta, nsample,\n"
            " mtype, order, dim, mic_orientation, hp_filter, stype, source_orientation);\n\n"
            "Input parameters:\n"
            " c                  : sound velocity in m/s.\n"
            " fs                 : sampling frequency in Hz.\n"
            " r                  : M x 3 array specifying the (x,y,z) coordinates of the\n"
            "                      receiver(s) in m.\n"
            " s                  : 1 x 3 vector specifying the (x,y,z) coordinates of the\n"
            "                      source in m.\n"
            " L                  : 1 x 3 vector specifying the room dimensions (x,y,z) in m.\n"
            " beta               : 1 x 6 vector specifying the reflection coefficients\n"
            "                      [beta_x1 beta_x2 beta_y1 beta_y2 beta_z1 beta_z2] or\n"
            "                      beta = reverberation time (T_60) in seconds.\n"
            " nsample            : number of samples to calculate, default is T_60*fs.\n"
            " mtype              : [omnidirectional, subcardioid, cardioid, hypercardioid,\n"
            "                      bidirectional], default is omnidirectional.\n"
            " order              : reflection order, default is -1, i.e. maximum order.\n"
            " dim                : room dimension (2 or 3), default is 3.\n"
            " mic_orientation    : direction in which the microphones are pointed, specified using\n"
            "                      azimuth and elevation angles (in radians), default is [0 0].\n"
            " hp_filter          : use 'false' to disable high-pass filter, the high-pass filter\n"
            "                      is enabled by default.\n\n"
            " stype              : [omnidirectional, subcardioid, cardioid, hypercardioid,\n"
            "                      bidirectional], default is omnidirectional.\n"
            " source_orientation : direction in which the source is pointed, specified using\n"
            "                      azimuth and elevation angles (in radians), default is [0 0].\n"
            "Output parameters:\n"
            " h                  : M x nsample matrix containing the calculated room impulse\n"
            "                      response(s).\n"
            " beta_hat           : In case a reverberation time is specified as an input parameter\n"
            "                      the corresponding reflection coefficient is returned.\n\n");
        return;
    }
    else
    {
        mexPrintf("Room Impulse Response Generator (Version 2.3.20241022) by Emanuel Habets\n"
            "Copyright (C) 2003-2025 E.A.P. Habets\n");
    }

    // Check for proper number of arguments
    if (nrhs < 6)
        mexErrMsgTxt("Error: There are at least six input parameters required.");
    if (nrhs > 14)
        mexErrMsgTxt("Error: Too many input arguments.");
    if (nlhs > 2)
        mexErrMsgTxt("Error: Too many output arguments.");

    // Check for proper arguments
    if (!(mxGetN(prhs[0])==1) || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("Invalid input arguments!");
    if (!(mxGetN(prhs[1])==1) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
        mexErrMsgTxt("Invalid input arguments!");
    if (!(mxGetN(prhs[2])==3) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
        mexErrMsgTxt("Invalid input arguments!");
    if (!(mxGetN(prhs[3])==3) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]))
        mexErrMsgTxt("Invalid input arguments!");
    if (!(mxGetN(prhs[4])==3) || !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]))
        mexErrMsgTxt("Invalid input arguments!");
    if (!(mxGetN(prhs[5])==6 || mxGetN(prhs[5])==1) || !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]))
        mexErrMsgTxt("Invalid input arguments!");

    // Load parameters
    double          c = mxGetScalar(prhs[0]);
    double          fs = mxGetScalar(prhs[1]);
    double*         rr = mxGetPr(prhs[2]);
    int             nMicrophones = (int) mxGetM(prhs[2]);
    double*         ss = mxGetPr(prhs[3]);
    double*         LL = mxGetPr(prhs[4]);
    double*         beta_input = mxGetPr(prhs[5]);
    double          beta[6];
    int             nSamples;
    char*           microphone_type;
    int             nOrder;
    int             nDimension;
    double          microphone_angle[2];
    int             isHighPassFilter;
    double          reverberation_time = 0;
    char*           source_type;
    double          source_angle[2];

    // Reflection coefficients or reverberation time?
    if (mxGetN(prhs[5])==1)
    {
        double V = LL[0]*LL[1]*LL[2];
        double S = 2*(LL[0]*LL[2]+LL[1]*LL[2]+LL[0]*LL[1]);
        reverberation_time = beta_input[0];
        if (reverberation_time != 0) {
            double alfa = 24*V*log(10.0)/(c*S*reverberation_time);
            if (alfa > 1)
                mexErrMsgTxt("Error: The reflection coefficients cannot be calculated using the current "
                             "room parameters, i.e. room size and reverberation time.\n           Please "
                             "specify the reflection coefficients or change the room parameters.");
            for (int i=0;i<6;i++)
                beta[i] = sqrt(1-alfa);
        }
        else
        {
            for (int i=0;i<6;i++)
                beta[i] = 0;
        }
    }
    else
    {
        for (int i=0;i<6;i++)
            beta[i] = beta_input[i];
    }

    // 3D Source orientation (optional)
    if (nrhs > 13 &&  mxIsEmpty(prhs[13]) == false)
    {
        double* source_orientation = mxGetPr(prhs[13]);
        if (mxGetN(prhs[13]) == 1)
        {
            source_angle[0] = source_orientation[0];
            source_angle[1] = 0;
        }
        else
        {
            source_angle[0] = source_orientation[0];
            source_angle[1] = source_orientation[1];
        }
    }
    else
    {
        source_angle[0] = 0;
        source_angle[1] = 0;
    }

    // Type of source (optional)
    if (nrhs > 12 &&  mxIsEmpty(prhs[12]) == false)
    {
        int return_value;
        source_type = new char[mxGetN(prhs[12])+1];        
        return_value = mxGetString(prhs[12], source_type, mxGetN(prhs[12])+1);
        if (return_value != 0)
        {
            mexErrMsgTxt("The input parameter stype is not a character array!");
        }
    }
    else
    {
        source_type = new char[1];
        source_type[0] = 'o';
    }
    
    // High-pass filter (optional)
    if (nrhs > 11 &&  mxIsEmpty(prhs[11]) == false)
    {
        isHighPassFilter = (int) mxGetScalar(prhs[11]);
    }
    else
    {
        isHighPassFilter = 1;
    }

    // 3D Microphone orientation (optional)
    if (nrhs > 10 &&  mxIsEmpty(prhs[10]) == false)
    {
        double* mic_orientation = mxGetPr(prhs[10]);
        if (mxGetN(prhs[10]) == 1)
        {
            microphone_angle[0] = mic_orientation[0];
            microphone_angle[1] = 0;
        }
        else
        {
            microphone_angle[0] = mic_orientation[0];
            microphone_angle[1] = mic_orientation[1];
        }
    }
    else
    {
        microphone_angle[0] = 0;
        microphone_angle[1] = 0;
    }

    // Room Dimension (optional)
    if (nrhs > 9 &&  mxIsEmpty(prhs[9]) == false)
    {
        nDimension = (int) mxGetScalar(prhs[9]);
        if (nDimension != 2 && nDimension != 3)
            mexErrMsgTxt("Invalid input arguments!");

        if (nDimension == 2)
        {
            beta[4] = 0;
            beta[5] = 0;
        }
    }
    else
    {
        nDimension = 3;
    }

    // Reflection order (optional)
    if (nrhs > 8 &&  mxIsEmpty(prhs[8]) == false)
    {
        nOrder = (int) mxGetScalar(prhs[8]);
        if (nOrder < -1)
            mexErrMsgTxt("Invalid input arguments!");
    }
    else
    {
        nOrder = -1;
    }

    // Type of microphone (optional)
    if (nrhs > 7 &&  mxIsEmpty(prhs[7]) == false)
    {
        int return_value;
        microphone_type = new char[mxGetN(prhs[7])+1];        
        return_value = mxGetString(prhs[7], microphone_type, mxGetN(prhs[7])+1);
        if (return_value != 0)
        {
            mexErrMsgTxt("The input parameter mtype is not a character array!");
        }
    }
    else
    {
        microphone_type = new char[1];
        microphone_type[0] = 'o';
    }

    // Number of samples (optional)
    if (nrhs > 6 &&  mxIsEmpty(prhs[6]) == false)
    {
        nSamples = (int) mxGetScalar(prhs[6]);
    }
    else
    {
        if (mxGetN(prhs[5])>1)
        {
            double V = LL[0]*LL[1]*LL[2];
            double alpha = ((1-pow(beta[0],2))+(1-pow(beta[1],2)))*LL[1]*LL[2] +
                ((1-pow(beta[2],2))+(1-pow(beta[3],2)))*LL[0]*LL[2] +
                ((1-pow(beta[4],2))+(1-pow(beta[5],2)))*LL[0]*LL[1];
            reverberation_time = 24*log(10.0)*V/(c*alpha);
            if (reverberation_time < 0.128)
                reverberation_time = 0.128;
        }
        nSamples = (int) (reverberation_time * fs);
    }

    // Create output vector
    plhs[0] = mxCreateDoubleMatrix(nMicrophones, nSamples, mxREAL);
    double* imp = mxGetPr(plhs[0]);

    // mexPrintf("microphone_type is %s; microphone_angle=(%f,%f)\n", microphone_type, microphone_angle[0] * 180/3.14159, microphone_angle[1] * 180/3.14159);
    // mexPrintf("source_type is %s; source_angle=(%f,%f)\n", source_type, source_angle[0] * 180/3.14159, source_angle[1] * 180/3.14159);
    // mexPrintf("c=%f \n fs=%f \n nMicrophones=%d \n nOrder=%d \n nDimension=%d \n", c, fs, nMicrophones, nOrder, nDimension);
    // mexPrintf("isHighPassFilter=%d \n reverberation_time=%f \n nSamples=%d \n", isHighPassFilter, reverberation_time, nSamples);
    computeRIR(imp, c, fs, rr, nMicrophones, nSamples, ss, LL, beta, microphone_type[0], nOrder, microphone_angle, isHighPassFilter, source_type[0], source_angle);

    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        double* beta_hat = mxGetPr(plhs[1]);
        if (reverberation_time != 0) {
            beta_hat[0] = beta[0];
        }
        else {
            beta_hat[0] = 0;
        }
    }

    delete[] microphone_type;
    delete[] source_type;
}
