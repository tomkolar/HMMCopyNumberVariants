/*
 * driver.cpp
 *
 *	This is the driver file for finding copy number variants
 *  using the maximal D-Segment algorithm.
 *
 *	Typical use:
 *		cnv cnvFile normalLength elevatedLength normalMean eleveatedMean
 *
 *  Created on: 3-16-13
 *      Author: tomkolar
 */
#include "DSegmentsFinder.h"
#include "HMMProbabilities.h"
#include <string>
#include <sstream>
#include <iostream>
using namespace std;

int main( int argc, char *argv[] ) {
/*
	// Check that file name, lengths and means were enetered as parameters
	if (argc < 5) {
			cout << "Invalid # of arguments\n";
			cout << "usage: cnv cnvFile normalLength elevatedLength normalMean eleveatedMean \n";
			return -1;
	}

	cout << "Starting\n";

	// Get Parameters
	string cnvFileName = argv[1];
	int normalLength = atoi(argv[2]);
	int elevatedLength = atoi(argv[3]);
	double normalMean = atof(argv[4]);
	double elevatedMean = atof(argv[5]);
*/
	// Set Parameters
	string cnvFileName = "C:/Users/kolart/Documents/Genome540/Assignment9/NA19238.chr20.counts";
	int normalLength = 1000000;
	int elevatedLength = 10000;
	double normalMean = 0.38;
	double elevatedMean = 0.57;

	// Create the DSegmentsFinder
	HMMProbabilities* probs = new HMMProbabilities(normalLength, elevatedLength, normalMean, elevatedMean);
	DSegmentsFinder* finder = new DSegmentsFinder(probs);
	cout << "D-Segments Finder Created.\n";

	// Find the d-segments
	finder->findDSegments(cnvFileName);
	cout << finder->results();

	cout << "Fred";
}
