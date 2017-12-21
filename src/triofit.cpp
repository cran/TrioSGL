
// For compilation testing: g++ -O3 -march=native -std=c++0x -ggdb3 -Wall -Wextra -pedantic -Wfatal-errors -fpic triofit.cpp -c -o triofit.o
//Remove -std=c++0x, if your compiler is too old

#include <cmath>
#include <algorithm>
#include <vector>
//#include <iostream>

using namespace std;

///Calculates the gradient of a given group, the result is stored in result
static void gradCalc(vector<double>& result, int groupStart, int groupLength, const vector<double>& minusSumX_i, const double* eta, const double* X, int nrow) {

    for (int i = 0; i < groupLength; ++i) {
        result[i] = minusSumX_i[i + groupStart];
    }

    vector<double> tempvector(groupLength);
    for (int i = 0; i < nrow; i += 4) {     //i is index of affected persons
        double tempdivisor[4] = {0, 0, 0, 0};
        for (int j = i; j < i + 4; ++j) {       //affected person and 3 pseudo controls
            for (int k = 0; k < 4; ++k) {
                tempdivisor[k] += exp(eta[j] - eta[i + k]);         //this formulation nullifies the impact of double-overflows in exp(some big value)
            }
        }
        for (int j = 0; j < groupLength; ++j) {
            for (int k = 0; k < 4; ++k) {
                result[j] += X[i + k + nrow * (j + groupStart)] / tempdivisor[k];
            }
        }
    }

    for (int i = 0; i < groupLength; ++i) {
        result[i] /= (nrow / 4);
    }

    return;
}


///Calculates the actual negative log-Likelihood
static double negLogLikelihoodCalc(const double* eta, int nrow) {

    double sum = 0;
    for (int i = 0; i < nrow; i += 4) {
        sum -= eta[i];
        double maxeta = *max_element(&eta[i], &eta[i + 4]);
        sum += maxeta;
        double triosum = 0;
        for (int j = i; j < i + 4; ++j) {
            triosum += exp(eta[j] - maxeta);
        }
        sum += log(triosum);
    }

    return sum / (nrow / 4);
}


///Uses gradient-based method on each group once to minimize SGL
static void Solver(const double* X, const vector<double>& minusSumX_i, int nrow, int ncol, int numGroup, double* beta, const int* rangeGroupInd, const int* groupLen, double lambda1, double lambda2, int innerIter, double thresh, double gamma, double* eta, int* betaIsZero, bool& groupChange, bool* isActive, bool* useGroup, double step, int reset) {
    double* theta = new double[ncol];
    double* etaNew = new double[nrow];

    const int maxGroupLen = *max_element(groupLen, groupLen + numGroup);
    vector<double> grad(maxGroupLen);       //reused for different grouplens

    for (int i = 0; i < numGroup; i++) {
        if (!useGroup[i]) continue;

        int startInd = rangeGroupInd[i];

        // Setting up null gradient calc to check if group is 0
        for (int k = 0; k < nrow; k++) {
            etaNew[k] = eta[k];
            for (int j = startInd; j < rangeGroupInd[i] + groupLen[i]; j++) {
                etaNew[k] -= X[k + nrow * j] * beta[j];
            }
        }

        // Calculating Null Gradient
        gradCalc(grad, rangeGroupInd[i], groupLen[i], minusSumX_i, etaNew, X, nrow);

        for (int j = 0; j < groupLen[i]; j++) {
            //Anwendung von S(...) wird berechnet, da später quadriert wird, muss sign(...) nicht beachtet werden und die Berechnung kann die Betragsfunktion anders berechnen
            if (grad[j] > lambda1) {
                grad[j] -= lambda1;
            } else if (grad[j] < -lambda1) {
                grad[j] += lambda1;
            } else { // ( lambda1 >= grad[j] && grad[j] >= -lambda1)
                grad[j] = 0;
            }
        }

        double zeroCheck = 0;
        for (int j = 0; j < groupLen[i]; j++) {
            zeroCheck += grad[j] * grad[j];
        }

        if (zeroCheck <= lambda2 * lambda2 * groupLen[i]) { //Or not? Formel passt mit S.6 unten, wenn sqrt(p_l) mit in den Koeff. eingerechnet ist
            if (betaIsZero[i] == 0) {
                for (int k = 0; k < nrow; k++) {
                    for (int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++) {
                        eta[k] -= X[k + nrow * j] * beta[j];
                    }
                }
                betaIsZero[i] = 1;
            }
            for (int j = 0; j < groupLen[i]; j++) {
                beta[j + rangeGroupInd[i]] = 0;
            }
        } else {
            if (!isActive[i]) {
                isActive[i] = true;
                groupChange = true;
            }

            for (int k = 0; k < ncol; k++) {
                theta[k] = beta[k];
            }

            betaIsZero[i] = 0;
            double* z = new double[groupLen[i]];
            double* U = new double[groupLen[i]];
            double* G = new double[groupLen[i]];
            double* betaNew = new double[ncol];

            int count = 0;
            double check = 1e100;

            while (count <= innerIter && check > thresh) {

                count++;

                gradCalc(grad, rangeGroupInd[i], groupLen[i], minusSumX_i, eta, X, nrow);

                double Lold = negLogLikelihoodCalc(eta, nrow);
                // Back-tracking

                while (true) {

                    for (int j = 0; j < groupLen[i]; j++) {
                        z[j] = beta[j + rangeGroupInd[i]] - step * grad[j];
                        if (z[j] <= -lambda1 * step) {
                            z[j] += lambda1 * step;
                        } else if (z[j] >= lambda1 * step) {
                            z[j] -= lambda1 * step;
                        } else {    //(z[j] < lambda1 * step && z[j] > -lambda1 * step)
                            z[j] = 0;
                        }
                    }

                    double norm = 0;
                    for (int j = 0; j < groupLen[i]; j++) {
                        norm += z[j] * z[j];
                    }
                    norm = sqrt(norm);

                    double uOp = 0;
                    if (norm != 0) {
                        uOp = (1 - lambda2 * sqrt(double(groupLen[i])) * step / norm);   //Or not?
                        if (uOp < 0) {
                            uOp = 0;
                        }
                    }

                    for (int j = 0; j < groupLen[i]; j++) {
                        U[j] = uOp * z[j];
                        G[j] = 1 / step * (beta[j + rangeGroupInd[i]] - U[j]);
                    }

                    for (int k = 0; k < nrow; k++) {
                        etaNew[k] = eta[k];
                        for(int j = 0; j < groupLen[i]; j++) {
                            etaNew[k] -= step * G[j] * X[k + nrow * (rangeGroupInd[i] + j)];
                        }
                    }

                    double Lnew = negLogLikelihoodCalc(etaNew, nrow);

                    double sqNormG = 0;
                    double iProd = 0;

                    for (int j = 0; j < groupLen[i]; j++) {
                        sqNormG += G[j] * G[j];
                        iProd += grad[j] * G[j];
                    }

                    //double diff = Lold - Lnew - step * iProd + step / 2 * sqNormG;
                    if (Lold - Lnew - step * iProd + step / 2 * sqNormG >= 0) {
                        break;
                    }

                    step *= gamma;
                }

                check = 0;

                for(int j = 0; j < groupLen[i]; j++) {
                    check += fabs(theta[j + rangeGroupInd[i]] - U[j]);
                    for(int k = 0; k < nrow; k++) {
                        eta[k] -= X[k + nrow * (j + rangeGroupInd[i])] * beta[j + rangeGroupInd[i]];
                    }
                    beta[j + rangeGroupInd[i]] = U[j] + (count % reset) / (count % reset + 3) * (U[j] - theta[j + rangeGroupInd[i]]);
                    theta[j + rangeGroupInd[i]] = U[j];

                    for(int k = 0; k < nrow; k++) {
                        eta[k] += X[k + nrow * (j + rangeGroupInd[i])] * beta[j + rangeGroupInd[i]];
                    }
                }
            }
            delete [] z;
            delete [] U;
            delete [] G;
            delete [] betaNew;
        }
    }
    delete [] etaNew;
    delete [] theta;
}


///Returns true, iff. |a| < |b|
static bool absLessCompare (double a, double b) {
    return fabs(a) < fabs(b);
}


///Encloses Solver in an outer loop with a max. iteration counter and threshold-based criterion
static void triofitOuterLoop (const double* X, int nrow, int ncol, int numGroup, const int* rangeGroupInd, const int* groupLen, double lambda1, double lambda2, double* beta, int innerIter, int outerIter, double thresh, double outerThresh, double* eta, double gamma, int* betaIsZero, double step, int reset) {

    bool* isActive = new bool[numGroup];
    bool* useGroup = new bool[numGroup];
    bool* tempIsActive = new bool[numGroup];

    for (int i = 0; i < numGroup; i++) {
        isActive[i] = false;
        useGroup[i] = true;
    }

    vector<double> minusSumX_i(ncol);

    for (int i = 0; i < nrow; i += 4) {
        for (int j = 0; j < ncol; ++j) {
            minusSumX_i[j] -= X[i + nrow * j];
        }
    }

    // outer most loop creating response etc...
    int outermostCounter = 0;   //outermostCounter is not reset in the outer loop, because outerIter shall be an overall maximal number of iterations
    double* outerOldBeta = new double[ncol];

    bool groupChange = true;
    while (groupChange) {
        groupChange = false;

        Solver(X, minusSumX_i, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, gamma, eta, betaIsZero, groupChange, isActive, useGroup, step, reset);

        while(outermostCounter < outerIter) {
            outermostCounter++;

            for (int i = 0; i < ncol; i++) {
                outerOldBeta[i] = beta[i];
            }

            for (int i = 0; i < numGroup; i++) {
                tempIsActive[i] = isActive[i];
            }

            Solver(X, minusSumX_i, nrow, ncol, numGroup, beta, rangeGroupInd, groupLen, lambda1, lambda2, innerIter, thresh, gamma, eta, betaIsZero, groupChange, isActive, tempIsActive, step, reset);

            double relativeErrorThreshold = fabs(*max_element(beta, beta + ncol, absLessCompare) * thresh * 0.1);       //Can be 0

            double outermostCheck = 0;
            for(int i = 0; i < ncol; i++) {
                if (fabs(beta[i]) > relativeErrorThreshold) {
                    outermostCheck += fabs((outerOldBeta[i] - beta[i]) / beta[i]);      //relative error
                } else {
                    outermostCheck += fabs(outerOldBeta[i] - beta[i]);                  //absolute error, if below relativeErrorThreshold
                }
            }
            if (outermostCheck < outerThresh) {
                break;
            }
        }
    }

    delete [] outerOldBeta;
    delete [] useGroup;
    delete [] isActive;
    delete [] tempIsActive;
}


extern "C" {
    ///Just forwarding the data to C++ and removing some of the pointers
    void triofit (double* X, int* nrow, int* ncol, int* numGroup, int* rangeGroupInd, int* groupLen, double* lambda1, double* lambda2, double* beta, int* innerIter, int* outerIter, double* thresh, double* outerThresh, double* eta, double* gamma, int* betaIsZero, double* step, int* reset) {
        triofitOuterLoop(X, *nrow, *ncol, *numGroup, rangeGroupInd, groupLen, *lambda1, *lambda2, beta, *innerIter, *outerIter, *thresh, *outerThresh, eta, *gamma, betaIsZero, *step, *reset);
    }
}

