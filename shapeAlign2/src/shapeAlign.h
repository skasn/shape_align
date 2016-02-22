#ifndef __SHAPEALIGN
#define __SHAPEALIGN

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <gsl/gsl_matrix.h>

using namespace std;

// Struct to store results of optimal alignment discovery
typedef struct alignData{
	double score;
	int shift;
	bool rev;
} alignData;

class shapeAlign{

public:
	// Constructor
	shapeAlign(const string& nameList, const vector<string> &files,
		const int &minS, const int &maxS, const bool &win, const int &wS,
		const int &wE, const bool &ign, const int &iS, const int &iE, const int &E);

	// Destructor
	~shapeAlign(void);

	// Output methods:
	void printCentroid(void);
	void printShifts(void);
	void printShiftMatrix(void);
	void printDistanceMatrix(void);
	void printRevMatrix(void);
	void printShiftedProfiles(void);

private:

	string nameFile;
	vector<string> shapeFiles;

	// Min and max shifts
	int shiftMin;
	int shiftMax;

	// Alignment window start and end
	bool window;
	int winStart;
	int winEnd;

	// Window to ignore
	bool ignore;
	int ignStart;
	int ignEnd;

	// Minimum number of valid bases in overlap
	int thresh;

	vector<string> names;
	vector<gsl_matrix*> matrices;

	// The following ars nSites x nSites (square) matrices that
	// summarize the pairwise distances (D), optimal shifts (S),
	// and reversal status (R). Reversal status refers to whether
	// the optimal alignment requires reversal of shape vectors
	// in a comparison.
	gsl_matrix* D;
	gsl_matrix* S;
	gsl_matrix* R;

	int nSites;
	int m;			// Number of shape parameters
	int cIdx;		// index of centroid sequence (in names vector)
	double cDist;	// Distance of centroid to other sites

	void scaleMatrix(gsl_matrix* M, double min, double max);
	void scaleMatrixZscore(gsl_matrix *M);
	double normFrobenius(const gsl_matrix* M);
	double cosineSim(const gsl_matrix* X, const gsl_matrix *Y);
	alignData getOptimalShift(gsl_matrix* A1, gsl_matrix* A2);
	void reverse(gsl_matrix *R);
	pair<int,double> getCentroid(void);
	void printMatrix(gsl_matrix* M, const string &outFile);
};

#endif
