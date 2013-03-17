/*
 * DSegmentFinder.h
 *
 *	This is the header file for the DSegmentFinder object. A DSegmentFinder
 *  can find the maximal D-Segments in a sequence.
 *
 *  Created on: 3-16-13
 *      Author: tomkolar
 */

#ifndef DSEGMENTFINDER_H
#define DSEGMENTFINDER_H
#include "HMMProbabilities.h"
#include <string>
#include <vector>
using namespace std;

class DSegmentsFinder
{
public:
	// Constuctors
	// ==============================================
	DSegmentsFinder();
	DSegmentsFinder(HMMProbabilities* probs);

	// Destructor
	// =============================================
	~DSegmentsFinder();

	// Public Attributes
	// =============================================
	HMMProbabilities* probabilities;

	// Public Methods
	// =============================================
	 
	// findDSegments(string cnvFileName)
	//  Purpose: 
	//		Finds the DSegments for the sequence
	void findDSegments(string cnvFileName);

	// string results()
	//  Purpose:
	//		Returns a string representing the results for finding the D-Segments
	//
	//		format:
	//			<result type="viterbi_iteration" iteration="<< iteration >>">
	//				<<probabilitiesResultsString>>
	//				<<thresholdResultsString>>
	//				<<segmentsResultsString>>
	//				<<readStartCountsAllResultsString>>
	//				<<readStartCountsDSegmentsResultsString>>
	//			</result>
	string results();

private:
	struct Segment {
		int start;
		int end;
		long double score;
	};

	vector<Segment> segments;
	int readStartCounts[4];
	int dSegmentReadStartCounts[4];
	double threshold;

	// string probabilitiesResultsString()
	//  Purpose:
	//		Returns a string representing the probabilites
	//
	//		format:
	//			<<statesResultsString>>
	//			<<initiationProbabilitesResultsString>>
	//			<<transmissionProbabilitesResultsString>>
	//			...
	//			<<emissionProbabilitesResultsString>>
	//			...
	string probabilitiesResultsString();

	// string thresholdResultsString()
	//  Purpose:
	//		Returns a string representing the threshold (S=-D)
	//
	//		format:
	//			<score_threshold>
	//				<<threshold>>,
	//			</score_threshold>
	string thresholdResultsString();

	// string segmentsResultsString()
	//  Purpose:
	//		Returns a string representing the segments
	//
	//		format:
	//			<result type="segments">
	//				(segment1start, segment1end, segement1Score),(segment2start, segment2end, segment2Score),...
	//			</result>
	string segmentsResultsString();

	// string readStartCountsAllResultsString()
	//  Purpose:
	//		Returns a string representing the read start counts
	//		for all positions in the sequence
	//
	//		format:
	//			<result type="read_start_counts_histogram" positions="all">
	//				<<#readStarts>>=<<readStartCount>>,...
	//			</result>
	string readStartCountsAllResultsString();

	// string readStartCountsDSegmentsResultsString()
	//  Purpose:
	//		Returns a string representing the read start counts
	//		for all positions in the D-Segments
	//
	//		format:
	//			<result type="read_start_counts_histogram" positions="state2">
	//				<<#readStarts>>=<<readStartCount>>,...
	//			</result>
	string readStartCountsDSegmentsResultsString();

};

#endif //DSEGMENTFINDER_H
