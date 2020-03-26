#ifndef _FILTER_H_
#define _FILTER_H_

#include <Rcpp.h>

RcppExport SEXP countCells (   //Inputs:
                        SEXP RFilter, // true if filter, false if thresholds
                        SEXP RPartitionsPerChannel, //no of partitions per channel
                        SEXP RThresholds,   //threshold points per channel; list of vectors
                        SEXP RmaxPopSize,
                        SEXP RthresChannels, // labels of threshold marker channels
                        SEXP RMfiData, // expr data from the flowFrame, only including mfi channels
                        SEXP RFrameExprData, //expr data from the flowFrame
                        SEXP RNumPops, //total number of populations
                        SEXP Rverbose //total number of populations
                        );


#endif // _FILTER_H_
