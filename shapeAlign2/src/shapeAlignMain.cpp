/*****************************************************************************
	simulateHORMain.cpp

	Siva Kasinathan
	Henikoff Lab
	Fred Hutchinson Cancer Research Center
	skasin@uw.edu

	Parameter checking macro/main borrowed from Aaron Quinlan's
	BEDTools suite.

******************************************************************************/
#include "shapeAlign.h"
#include <algorithm>
#include <cstring>
#include <string>
#include <sstream>

using namespace std;

// define our program name
#define PROGRAM_NAME "shapeAlign"
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, std::min(actualLen, paramLen))== 0) && (actualLen == paramLen)


// function declarations
void shapeAlign_help(void);

int main(int argc, char* argv[]) {

    bool showHelp = false;
    bool haveShapeFiles = false;
    bool haveNameFile = false;
    bool window = false;
    bool ignore = false;

	string namesFile;
	vector<string> shapeFiles;

    // Default max and min shifts
	int max = 25;
	int min = -25;
	int start = -1;
	int end = -1;
    int istart = -1;
    int iend = -1;
    int entr = 10;

    // check to see if we should print out some help
    if(argc <= 1) showHelp = true;
    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }
    if(showHelp) shapeAlign_help();

     for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-n",2, parameterLength)) {
            if ((i+1) < argc) {
                haveNameFile = true;
				namesFile = argv[i+1];
                i++;
            }
        }
        else if(PARAMETER_CHECK("-f",2, parameterLength)) {
            haveShapeFiles = true;
            while ((i+1) < argc) {
				shapeFiles.push_back(argv[i+1]);
                i++;
            }
        }
        else if(PARAMETER_CHECK("-min",4, parameterLength)) {
            if ((i+1) < argc) {
            	stringstream ss(argv[i+1]);
            	ss >> min;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-max",4, parameterLength)) {
            if ((i+1) < argc) {
            	stringstream ss(argv[i+1]);
            	ss >> max;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-start",6, parameterLength)) {
            if ((i+1) < argc) {
            	stringstream ss(argv[i+1]);
            	ss >> start;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-end",4, parameterLength)) {
            if ((i+1) < argc) {
            	stringstream ss(argv[i+1]);
            	ss >> end;
                i++;
            }
        }
         else if(PARAMETER_CHECK("-istart",7, parameterLength)) {
            if ((i+1) < argc) {
                stringstream ss(argv[i+1]);
                ss >> istart;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-iend",5, parameterLength)) {
            if ((i+1) < argc) {
                stringstream ss(argv[i+1]);
                ss >> iend;
                i++;
            }
        }
        else if(PARAMETER_CHECK("-entr",5, parameterLength)) {
            if ((i+1) < argc) {
                stringstream ss(argv[i+1]);
                ss >> entr;
                i++;
            }
        }
        else {
            cerr << "*****ERROR: Unrecognized parameter: "
                 << argv[i]
                 << " *****"
                 << endl << endl;
            showHelp = true;
        }
    }


	if (!haveShapeFiles || !haveNameFile) {
        showHelp = true;
    }

    if (start != -1 && end != -1){
        if (start <= end) {window = true; }
        else {
            cerr << "Alignment window start exceeds alignment window end. Ignoring alignment window." << endl;
        }
    } else if ((start != -1 || end != -1) && (start == -1 || end == -1)){
        cerr << "Must define BOTH alignment start and end if -start/-end" << endl;
        cerr << "flags are used." << endl << endl;
        showHelp=true;
    }

    if (istart != -1 && iend != -1){
        if (istart <= iend) {ignore = true; }
        else {
            cerr << "Ignore position start exceeds ignore position end." << endl << "Considering all positions in scoring." << endl;
        }
    } else if ((istart != -1 || iend != -1) && (istart == -1 || iend == -1)){
        cerr << "Must define BOTH ignore position start and end if -istart/-iend" << endl;
        cerr << "flags are used." << endl << endl;
        showHelp=true;
    }

    if (!showHelp) {
        shapeAlign *align = new shapeAlign(namesFile,shapeFiles,min,max,window,start,end,ignore,istart,iend,entr);
        delete align;
    }
    else {
    	shapeAlign_help();
    }


    return 0;
}

void shapeAlign_help(void) {

    cerr << "\nSequence-independent alignment of DNA shape" << endl << endl;

    cerr << "Usage:   " << PROGRAM_NAME
         << " [OPTIONS] -n <site name file> -f <shape file 1> ... <shape file n> "
         << endl << endl;
    cerr << " *Note: shape file list must be the last entered command line argument."
         << endl << endl;
    cerr << " -n      Single-column file containing names/accessions of sites." << endl;
    cerr << " -f      Space-delimited list of shape files (multi-column, tab-" << endl;
    cerr << "         delimited file). By default, the first three columns and" << endl;
    cerr << "         the last two columns are ignored. The ordering of sites" << endl;
    cerr << "         in these files is assumed to be the same as in the site" << endl;
    cerr << "         names file." << endl << endl;
    cerr << " Options: " << endl;
    cerr << " -entr   Minimum number of valid bases in overlap. (Default: 10)." << endl;
    cerr << " -min    Minimum shift. (Default: -25 bp)." << endl;
    cerr << " -max    Maximum shift. (Default: 25 bp)." << endl;
    cerr << " -start  Start of alignment window defined with respect to the" << endl;
    cerr << "         start of the shape window, which is defined as zero." << endl;
    cerr << "         (Note: If -start is used, -end must be defined)." << endl;
    cerr << " -end    End of alignment window defined with respect to the " << endl;
    cerr << "         start of the shape window, which is defined as zero." << endl;
    cerr << "         (Note: If -end is used, -start must be defined)." << endl;
    cerr << " -istart  Start of positions to ignore defined with respect to" << endl;
    cerr << "         start of the shape window, which is defined as zero." << endl;
    cerr << "         (Note: If -istart is used, -iend must be defined)." << endl;
    cerr << " -iend    End of positions to ignore defined with respect to" << endl;
    cerr << "         start of the shape window, which is defined as zero." << endl;
    cerr << "         (Note: If -iend is used, -istart must be defined)." << endl << endl;
    exit(1);

}
