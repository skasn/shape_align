#include "shapeAlign.h"
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <ostream>
#include <sstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <omp.h>

// Constructor
shapeAlign::shapeAlign(const string& nameList, const vector<string> &files,
	const int &minS, const int &maxS, const bool &win, const int &wS,
	const int &wE, const bool &ign, const int &iS, const int &iE, const int &E):
	nameFile(nameList), shapeFiles(files), shiftMin(minS), shiftMax(maxS),
	window(win), winStart(wS), winEnd(wE), ignore(ign), ignStart(iS), ignEnd(iE),
	thresh(E)
{

	// Get number of shape parameters -- one file for each parameter,
	// so this is effectively the number of files
	m = files.size();

	// Get site names by reading the single-column file containing
	// names of sites
	ifstream file(nameFile.c_str());
	string line;

	while(getline(file,line))
		names.push_back(line);
	file.close();

	// Get number of sites
	nSites = names.size();

	cerr << "Read " << nSites << " site names." << endl;

	// Initialize matrices for tracking pairwise comparison info
	D = gsl_matrix_calloc(nSites,nSites);
	S = gsl_matrix_calloc(nSites,nSites);
	R = gsl_matrix_calloc(nSites,nSites);

	cerr << "Reading shape files." << endl;

	matrices.resize(nSites);	// Initialize an empty list of matrix references

	// create matrices containing the shape information. Add data
	// to these matrices on the fly.
	for (size_t f = 0; f < files.size(); f++){

		cerr << "\t" << files[f] << endl;

		ifstream shapeFile(files[f].c_str());
		int idx = 0;	// line/site counter
		while(getline(shapeFile,line)){		// Each line in the shape files represent a single site
			stringstream linestream(line);
			string s;
			vector <string> temp;
			while(linestream >> s)	// Split on spaces and store data from each position in site
				temp.push_back(s);

			int n = temp.size();	// Get number of positions

			// there are five columns that need to be trimmed off:
			// The first three columns (identifier and NAs) and the
			// last two columns (NAs)

			// Initialize the matrix if the matrix has not previously been
			// initialized
			if (f == 0)	matrices[idx] = gsl_matrix_alloc(m,n-5);

			for (size_t i = 0; i < matrices[idx]->size2; i++){
				double d;
				stringstream stod(temp[i+3]);
				stod >> d;
				gsl_matrix_set(matrices[idx],f,i,d);
				}
			// Increment the line counter
			idx++;
		}
	}


	cerr << "\tDone reading shape files." << endl;

	// Scale each matrix such that values to go from 1->2
	cerr << "Scaling matrices." << endl;

	for (size_t i = 0; i < nSites; i++)
		scaleMatrixZscore(matrices[i]);

	cerr << "\tDone scaling matrices." << endl;

	// Loop over the sites and compute all pairwise distances -- note that
	// distances are symmetric: D[a,b] = D[b,a]. But, the shifts computed
	// are not symmetric: S[a,b] = -S[b,a].
	for (size_t i = 0; i < nSites; i++){

		if ((i+1) % 100 == 0)
			cerr << "\tProcessing " << i+1 << " of " << nSites << endl;

		// Parallelize this portion: data races shouldn't be a concern
		// since no threads should be writing to the same block of
		// memory

		#pragma omp parallel
		{
			#pragma omp master
			if (i==0)
				cerr << "Beginning all-by-all distance calculation using "
					 << omp_get_num_threads() << " threads." << endl;


			#pragma omp for
			for (size_t j = i; j < nSites; j++){

				// Get optimal shift and distances for the simple comparison
				alignData results = getOptimalShift(matrices[i],matrices[j]);

				// Get the matrix representing the reverse of the matrices[j]
				gsl_matrix* rev = gsl_matrix_alloc(matrices[j]->size1,matrices[j]->size2);
				gsl_matrix_memcpy(rev,matrices[j]);
				reverse(rev);

				// Get the optimal shift and distance for the reverse matrix
				alignData resultsRev = getOptimalShift(matrices[i],rev);

				if (results.score >= resultsRev.score){
					results.rev = 0;
				} else {
					results.score = resultsRev.score;
					results.shift = resultsRev.shift;
					results.rev = 1;
				}

				// Store the data in the matrices used for tracking
				// pairwise comparisons
				gsl_matrix_set(D,i,j,results.score);
				gsl_matrix_set(S,i,j,results.shift);
				gsl_matrix_set(R,i,j,results.rev);

				gsl_matrix_set(D,j,i,results.score);
				gsl_matrix_set(S,j,i,-1*results.shift);
				gsl_matrix_set(R,j,i,results.rev);

				// Clean up -- free memory associated with rev
				gsl_matrix_free(rev);
			}
		}
	}

	cerr << "\tDone with distance calculation." << endl;

	cerr << "Finding centroid." << endl;
	pair<int,double> C = getCentroid();
	cIdx = C.first;		// Index (w.r.t. names vector) of centroid
	cDist = C.second;	// Distance of centroid to other sequences
	cerr << "\tCentroid: Site \"" << names[C.first] << "\"" << endl;
	cerr << "\tDistance: "  << C.second << endl;
	printCentroid();
	printShifts();

// 	cerr << "Printing matrices to files." << endl;
// 	printShiftMatrix();
	// printDistanceMatrix();
// 	printRevMatrix();
// 	cerr << "\tDone." << endl;

	cerr << "Printing aligned data." << endl;
	printShiftedProfiles();
	cerr << "\tDone." << endl;

	cerr << "Job successfully completed." << endl;

}

// Destructor
shapeAlign::~shapeAlign(void){
	// Free GSL matrices
	for (size_t i = 0; i < matrices.size(); i++)
		gsl_matrix_free(matrices[i]);
}

// Reverse the columns of a matrix
void shapeAlign::reverse(gsl_matrix* rev){
	for (size_t i = 0; i < rev->size2/2; i++)
		gsl_matrix_swap_columns(rev,i,rev->size2-i-1);
	return;
}

// For a number of horizontal shifts, find the optimal shift
alignData shapeAlign::getOptimalShift(gsl_matrix *A1, gsl_matrix *A2){
	alignData results;			// Container for holding optimal shift results
	results.score = -10000;		// Initialize results.score to a large negative number
	results.shift = 0;
	double score;

	// Perform a bunch of shifts, keeping track of the optimal score
	for (int s = shiftMin; s <= shiftMax; s++){
		// Initialize a matrix delta to contain the element-wise differences
		// in the overlapping region
		int delSize = A2->size2-abs(s);

		// Get the values that fall within the overlapping region as
		// submatrices
		gsl_matrix_view subA1v = (s <= 0 ? gsl_matrix_submatrix(A1,0,0,m,A1->size2+s)
			: gsl_matrix_submatrix(A1,0,s,m,A1->size2-s) );
		gsl_matrix_view subA2v = ( s <= 0 ? gsl_matrix_submatrix(A2,0,abs(s),m,A2->size2+s)
			: gsl_matrix_submatrix(A2,0,0,m,A2->size2-s) );


		gsl_matrix *S1 = gsl_matrix_alloc(m,delSize);
		gsl_matrix_memcpy(S1,&subA1v.matrix);
		gsl_matrix *S2 = gsl_matrix_alloc(m,delSize);
		gsl_matrix_memcpy(S2,&subA2v.matrix);

		for (size_t j = 0; j < delSize; j++){
			for (size_t i = 0; i < m; i++){
				if (window && (j <= winStart || j >= winEnd))
					s <=0 ? gsl_matrix_set(S1,i,j,0) : gsl_matrix_set(S2,i,j,0);
				if (ignore && (j >= ignStart && j <= ignEnd))
					s <=0 ? gsl_matrix_set(S1,i,j,0) : gsl_matrix_set(S2,i,j,0);
			}
		}

		score = cosineSim(S1,S2);

		if (score > results.score ){
			results.score = score;
			results.shift = s;
		}

		gsl_matrix_free(S1);
		gsl_matrix_free(S2);

	}
	return results;
}

// Calculate the 'normalized' Frobenius norm (i.e., divide
// by the square root of the number of elements in matrix).
// Instead of doing element-wise computation, use the GSL
// BLAS interface to compute the L2 norm of the columns
// and get the Frobenius norm as the sum of the squares of
// L2 norms of the column vectors
double shapeAlign::normFrobenius(const gsl_matrix *M)
{
	double L2n, Fn = 0;
	double val = 0;

	for (size_t i = 0; i  < M->size2; i++){
		if (! gsl_isnan(gsl_matrix_get(M,0,i))){
		    gsl_vector_const_view column = gsl_matrix_const_column(M,i);
		    L2n = gsl_blas_dnrm2(&column.vector);
	  	  Fn += L2n*L2n;
	  	  val+=1.0;
	  	 }
	}
	if (val >= 4*thresh){
		return sqrt(Fn/val);
	}
	return GSL_NAN;
}

double shapeAlign::cosineSim(const gsl_matrix* X, const gsl_matrix *Y)
{
	double cosSim = 0;
	double denomX = 0, denomY = 0;
	for (size_t i = 0; i < X->size1; i++){
		double dot = 0;
		// Get vector views of rows
		gsl_vector_const_view Xrow = gsl_matrix_const_row(X,i);
		gsl_vector_const_view Yrow = gsl_matrix_const_row(Y,i);

		// Compute dot product
		gsl_blas_ddot(&Xrow.vector,&Yrow.vector,&dot);
		cosSim +=dot;

		denomX += pow(gsl_blas_dnrm2(&Xrow.vector),2.0);
		denomY += pow(gsl_blas_dnrm2(&Yrow.vector),2.0);
	}
	cosSim = cosSim / (sqrt(denomX * denomY));
	return cosSim;
}

// Scale each row in a matrix to have values in a given range [a,b]
// f(x) = a + (b-a)(x-min)/(max-min)
void shapeAlign::scaleMatrix(gsl_matrix *M, double min, double max){
	for (size_t i = 0; i < M->size1; i++){
		double rmin,rmax;	// row max and min
		gsl_vector_view row = gsl_matrix_row(M,i);
		gsl_vector_minmax(&row.vector,&rmin,&rmax);
		gsl_vector_add_constant(&row.vector,-1*rmin);
		gsl_vector_scale(&row.vector,(max-min)/(rmax-rmin));
		gsl_vector_add_constant(&row.vector,min);
	}
	return;
}

void shapeAlign::scaleMatrixZscore(gsl_matrix *M){
	for (size_t i = 0; i < M->size1; i++){
		gsl_vector_view row = gsl_matrix_row(M,i);
		double mu = gsl_stats_mean(row.vector.data, row.vector.stride, row.vector.size);
		double sigma = gsl_stats_sd_m(row.vector.data, row.vector.stride, row.vector.size, mu);

		gsl_vector_add_constant(&row.vector,-mu);
		gsl_vector_scale(&row.vector,1.0/sigma);
	}
	return;
}

// Get the centroid, i.e., the site that has the lowest
// distance to other sites. This is done by computing the
// L2 norm (using the GLS BLAS interface) of the columns of
// the distance matrix. This is effectively the distance of
// a given sequence to all other sequences. The minimum
// value of the L2 norms gives corresponds to the centroid
pair<int,double> shapeAlign::getCentroid(void){
	gsl_vector * dists = gsl_vector_alloc(nSites);

	gsl_matrix_add_constant(D,-1.0);

	for (size_t i = 0; i < nSites; i++){
		gsl_vector_view column = gsl_matrix_column(D,i);
		gsl_vector_set(dists,i,gsl_blas_dnrm2(&column.vector));
	}

	gsl_matrix_add_constant(D,1.0);

	int minIdx = gsl_vector_min_index(dists);
	// int minIdx = gsl_vector_max_index(dists);
	double dist = gsl_vector_get(dists,minIdx);
	gsl_vector_free(dists);
	return make_pair(minIdx,dist);
}

void shapeAlign::printCentroid(void){
	int lastidx = nameFile.find_last_of(".");
	string centroidOutFile = nameFile.substr(0,lastidx) + ".centroid.txt";
	ofstream centroidInfo (centroidOutFile.c_str());

	if (centroidInfo.is_open()){
		centroidInfo << "Centroid: " << names[cIdx] << endl;
		centroidInfo << "Distance: " << cDist << endl;
	} else {
		cerr << "Error: Could not open output file ("
			 << centroidOutFile << "). Exiting." << endl;
		exit(1);	}
	centroidInfo.close();

	return;
}

void shapeAlign::printShifts(void){
	int mult;

	int lastidx = nameFile.find_last_of(".");
	string shiftOutFile = nameFile.substr(0,lastidx) + ".shift.stats.txt";
	ofstream shiftData (shiftOutFile.c_str());

	if (shiftData.is_open()){
		// Loop over all sites and output data from shape alignment
		for (size_t i = 0; i < nSites; i++){

			// Adjust the output for reversed sequences -- note that this
			// adjustment is only necessary in the summary output but not in
			// the actual shifted output!
			gsl_matrix_get(R,i,cIdx) ? mult = 1 : mult = -1;

			shiftData << names[i] << "\t" << gsl_matrix_get(D,i,cIdx)
		    	 << "\t" << mult*gsl_matrix_get(S,i,cIdx) << "\t"
		      	 << gsl_matrix_get(R,i,cIdx) << "\t"
		      	 << endl;
		}
	} else {
		cerr << "Error: Could not open output file ("
			 << shiftOutFile << "). Exiting." << endl;
		exit(1);
	}

	shiftData.close();

	return;
}

// Generic method for printing a matrix provided an output
// file name
void shapeAlign::printMatrix(gsl_matrix* M, const string &outFile){
	ofstream matrixOut(outFile.c_str());

	if (matrixOut.is_open()){
		for (size_t i = 0; i < M->size1; i++){
			for (size_t j = 0; j < M->size2; j++){
				if (j > 0) matrixOut << "\t";
				matrixOut << gsl_matrix_get(M,i,j);
			}
			matrixOut << endl;
		}
	} else {
		cerr << "Error: Could not open output file ("
			 << outFile << "). Exiting." << endl;
		exit(1);
	}

	matrixOut.close();
	return;
}

void shapeAlign::printShiftMatrix(void){
	int lastidx = nameFile.find_last_of(".");
	string shiftMatrixFile = nameFile.substr(0,lastidx) + ".shift.matrix.txt";

	printMatrix(S,shiftMatrixFile);

	return;
}

void shapeAlign::printDistanceMatrix(void){
	int lastidx = nameFile.find_last_of(".");
	string distMatrixFile = nameFile.substr(0,lastidx) + ".dist.matrix.txt";

	printMatrix(D,distMatrixFile);
	return;
}

void shapeAlign::printRevMatrix(void){
	int lastidx = nameFile.find_last_of(".");
	string revMatrixFile = nameFile.substr(0,lastidx) + ".rev.matrix.txt";

	printMatrix(R,revMatrixFile);
	return;
}

// Print shifted shape parameters. Given that the data read in at the
// outset were scaled, data will need to be read in again and appropriately
// shifted. This can be done on-the-fly given that the shifts have been
// computed and stored in the matrix S.
void shapeAlign::printShiftedProfiles(void){

	// Loop over all shape files
	for (size_t i = 0; i < m; i++){

		// Generate an output file name
		int lastidx = shapeFiles[i].find_last_of(".");
		string shapeOutFile = shapeFiles[i].substr(0,lastidx) + ".aligned" +
			shapeFiles[i].substr(lastidx,shapeFiles[i].length());

		ofstream shapeOut(shapeOutFile.c_str());

		if (shapeOut.is_open()){

			// Open the original shape data file
			ifstream shapeFile(shapeFiles[i].c_str());
			int idx = 0;
			string line;

			// Loop through all sites in the original shape data file
			while(getline(shapeFile,line)){
				stringstream linestream(line);
				string s;
				vector <string> temp;
				string name;
				int ctr = 0;	// Column counter
				while(linestream >> s){
					// Clip the first three columns
					if (ctr==0) name = s;
					if (ctr >2) temp.push_back(s);
					ctr++;
				}

				// Get rid of the last two columns
				temp.pop_back();
				temp.pop_back();

				// Get the optimal shift stored in S
				int shift = gsl_matrix_get(S,idx,cIdx);

				// If the shape data need to be reversed w.r.t. centroid,
				// do this now
				if (gsl_matrix_get(R,idx,cIdx))
					std::reverse(temp.begin(),temp.end());

				// Output data -- again, the first three columns and the
				// last two columns need to be trimmed because these
				// contain non-numeric information. There are two cases
				// depending on whether the columns were reversed or not

				shapeOut << name << "\t";

				for (size_t j = 0; j < temp.size(); j++){
					if (j > 0) shapeOut << "\t";

					if (j + shift >= 0 && j + shift < temp.size()){
						shapeOut << temp[j+shift];
					} else {
						shapeOut << "NA";
					}
				}
				shapeOut << endl;

				idx++;
			}


			shapeFile.close();

		}  else {
			cerr << "Error: Could not open output file ("
				 << shapeOutFile << "). Exiting." << endl;
			exit(1);
		}
		shapeOut.close();
	}

	return;
}

