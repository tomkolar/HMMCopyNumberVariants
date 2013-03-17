/*
 * DSegmentFinder.cpp
 *
 *	This is the cpp file for the DSegmentFinder object. A DSegmentFinder
 *  can find the maximal D-Segments in a sequence.
 *
 *  Created on: 3-16-13
 *      Author: tomkolar
 */
#include "DSegmentsFinder.h"
#include "StringUtilities.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

DSegmentsFinder::DSegmentsFinder() {
}

DSegmentsFinder::DSegmentsFinder(HMMProbabilities* probs) {
	// Initiailze counts
	for (int i = 0; i < 4; i++) {
		readStartCounts[i] = 0;
		dSegmentReadStartCounts[i] = 0;
	}

	// Initialize the probabailities and threshold
	probabilities = probs;
/*	threshold =
		log(
			(probs->transitionProbability(1,1) * probs->transitionProbability(2,2))
			/
			(probs->transitionProbability(1,2) * probs->transitionProbability(2,1))
		)
		/ log(2);
*/
		long double  sameSegProb =
			probs->logTransitionProbability(1,1) + probs->logTransitionProbability(2,2);
		long double  switchSegProb =
			probs->logTransitionProbability(1,2) + probs->logTransitionProbability(2,1);
		threshold = (sameSegProb - switchSegProb) / log(2);
}

DSegmentsFinder::~DSegmentsFinder() {
}

// findDSegments(string cnvFileName)
//  Purpose: 
//		Finds the DSegments for the sequence
void DSegmentsFinder::findDSegments(string cnvFileName) {

	ifstream inputFile(cnvFileName);
	string line;

	// Current segment read start counts vector
	int currentSegmentReadStartCounts[4];
	for (int i = 0; i < 4; i++)
		currentSegmentReadStartCounts[i] = 0;

	long double cum = 0;
	long double max = 0;
	int start = 1;
	int end = 1;
	int position, readStarts;
	while(getline(inputFile, line)) {

		//  Split the line into tokens to get the positon and readStarts
		vector<string> tokens;
		StringUtilities::split(line, '\t', tokens);
		position = atoi(tokens[1].c_str());
		readStarts = atoi(tokens[2].c_str());

		// Set read starts to 3 if greater than 3
		if (readStarts > 3)
			readStarts = 3;
	
		// Increment read start counts and temp counts;
		readStartCounts[readStarts]++;
		currentSegmentReadStartCounts[readStarts]++;

		// Add the score to the cumulative score
		cum += probabilities->dSegmentScore(readStarts);

		// Keep track of maximum score to this point
		if (cum >= max) {
			max = cum;
			end = position;
		}

		// Check if over threshold
		if (cum <= 0 || cum <= max - threshold) {
			if (max >= threshold) {
				// Create segment and add to segments collection
				Segment segment;
				segment.start = start;
				segment.end = end;
				segment.score = max;
				segments.push_back(segment); 

				// Add current segment counts to d-segment counts
				for (int i = 0; i < 4; i++)
					dSegmentReadStartCounts[i]+=currentSegmentReadStartCounts[i];
			}

			// Reset values
			cum = 0;
			max = 0;
			start = position + 1;
			end = position + 1;
			for (int i = 0; i < 4; i++)
				currentSegmentReadStartCounts[i] = 0;
		}
	}

	// Check if last segment is a D-Segment
	if (max >= threshold) {
		// Create segment and add to segments collection
		Segment segment;
		segment.start = start;
		segment.end = end;
		segment.score = max;
		segments.push_back(segment); 

		// Add current segment counts to d-segment counts
		for (int i = 0; i < 4; i++)
			dSegmentReadStartCounts[i]+=currentSegmentReadStartCounts[i];
	}

}

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
string DSegmentsFinder::results() {
	stringstream ss;

	// Header
	ss << "  <results>\n";

	// Results
	ss 
		<< probabilitiesResultsString()
		<< thresholdResultsString()
		<< segmentsResultsString()
		<< readStartCountsAllResultsString()
		<< readStartCountsDSegmentsResultsString();

	// Footer
	ss << "  </results>\n";

	return ss.str();
}

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
string DSegmentsFinder::probabilitiesResultsString() {
	return probabilities->probabilitiesResultsString();
}

// string thresholdResultsString()
//  Purpose:
//		Returns a string representing the threshold (S=-D)
//
//		format:
//			<score_threshold>
//				<<threshold>>,
//			</score_threshold>
string DSegmentsFinder::thresholdResultsString() {
	stringstream ss;
	ss.scientific;

	// Header
	ss << "      <score_threshold>";

	// Results
	ss 
		<< threshold;

	// Footer
	ss << "      </score_threshold>\n";

	return ss.str();
}


// string segmentsResultsString()
//  Purpose:
//		Returns a string representing the segments
//
//		format:
//			<result type="segments">
//				(segment1start, segment1end, segement1Score),(segment2start, segment2end, segment2Score),...
//			</result>
string DSegmentsFinder::segmentsResultsString() {
	stringstream ss;
	
	int counter = 0;
	for (int i = 0; i < segments.size();  i++) {
		Segment segment = segments[i];

		// Round score to one decimal place
		double score = segment.score;
		score = (score * 10) + .05;
		score = floor(score);
		score = score / 10;

		ss 
			<< "("
			<< segment.start
			<< ","
			<< segment.end
			<< ","
			<< score
			<< ")";
		if (i > 0)
			ss << ",";

		counter++;
		if (counter % 5 == 0)
			ss << "\n";
	}

	return StringUtilities::xmlResult("segment_list", ss.str());
}

// string readStartCountsAllResultsString()
//  Purpose:
//		Returns a string representing the read start counts
//		for all positions in the sequence
//
//		format:
//			<result type="read_start_counts_histogram" positions="all">
//				<<#readStarts>>=<<readStartCount>>,...
//			</result>
string DSegmentsFinder::readStartCountsAllResultsString() {
	stringstream ss;

	// Header
	ss << "    <result type=\"read_start_counts_histogram\" positions=\"all\">\n";

	// Results
	ss	<< "      ";
	for (int i = 0; i < 4; i++) {
		ss 
			<< i
			<< "="
			<< readStartCounts[i];

		if (i < 3)
			ss << ", ";
	}
	ss	<< "\n";

	// Footer
	ss << "    </result>\n";

	return ss.str();
}

// string readStartCountsDSegmentsResultsString()
//  Purpose:
//		Returns a string representing the read start counts
//		for all positions in the D-Segments
//
//		format:
//			<result type="read_start_counts_histogram" positions="state2">
//				<<#readStarts>>=<<readStartCount>>,...
//			</result>
string DSegmentsFinder::readStartCountsDSegmentsResultsString() {
	stringstream ss;

	// Header
	ss << "    <result type=\"read_start_counts_histogram\" positions=\"state2\">\n";

	// Results
	ss	<< "      ";
	for (int i = 0; i < 4; i++) {
		ss 
			<< i
			<< "="
			<< dSegmentReadStartCounts[i];

		if (i < 3)
			ss << ", ";
	}
	ss	<< "\n";

	// Footer
	ss << "    </result>\n";

	return ss.str();

}
