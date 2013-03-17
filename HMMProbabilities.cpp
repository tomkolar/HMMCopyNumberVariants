/*
 * HMMProbablities.cpp
 *
 *	This is the cpp file for the HMMProbabilties object. HMMProbabilities
 *  is a collection of all probabilties needed for a hidden markov model. This
 *  includes initiation, emission and transition probabilties.  There are 
 *	convenience methods for setting and retriving probabilties as well as
 *  the log value of each probabilty.
 *
 *  Created on: 2-15-13
 *      Author: tomkolar
 */
#include "HMMProbabilities.h"
#include "StringUtilities.h"
#include <cmath>
#include <sstream>
#include <string>
#include <limits>
#include <iostream>
#include <fstream>

// Constuctors
// ==============================================
HMMProbabilities::HMMProbabilities() {
}

HMMProbabilities::HMMProbabilities(int numOfStates) {
	numStates = numOfStates;
	createEmissionResidueMap();

	// Initialize all probabilities to zero
	for (int i = 0; i < numStates; i++) {
		setInitiationProbability(i, 0);
		for (int j = 0; j < numStates; j++) {
			setTransitionProbability(i, j, 0);
		}
		for (pair<string, int> mapPair : emissionResidueMap) {
			string& residue = mapPair.first;
			setEmissionProbability(i, residue, 0);
		}
	}
}

HMMProbabilities::HMMProbabilities(int normalLength, int elevatedLength, double normalMean, double elevatedMean) {
	numStates = 3;
	createEmissionResidueMap();
	setTransitionProbability(1, 1, 1 - ((double) 1/ (double) normalLength));
	setTransitionProbability(1, 2,  (double) 1/(double) normalLength);
	setTransitionProbability(2, 1,  (double) 1/(double) elevatedLength);
	setTransitionProbability(2, 2, 1 - ((double) 1/ (double) elevatedLength));

	populateEmissionProbabilities(1, normalMean);
	populateEmissionProbabilities(2, elevatedMean);
}


// Destructor
// =============================================
HMMProbabilities::~HMMProbabilities() {
}

// Public Methods
// =============================================

// double emissionProbability(int state, char residue)
//  Purpose: 
//		Returns the emission probability for the state and residue
long double HMMProbabilities::emissionProbability(int state, string residue) {
	return emissionProbabilities.at(state).at(getEmissionResidueIndex(residue));
}

// double initiationProbability(int state)
//  Purpose: 
//		Returns the initiation probability for the state
long double HMMProbabilities::initiationProbability(int state) {
	return initiationProbabilities[state];
}
	
// double transitionProbability(int beginState, int endState)
//  Purpose: 
//		Returns the transition probability for transition from beginState
//		to endState
long double HMMProbabilities::transitionProbability(int beginState, int endState) {
	return transitionProbabilities[beginState][endState];
}

// double logEmissionProbability(int state, char residue)
//  Purpose: 
//		Returns the log of the emission probability for the state and residue
long double HMMProbabilities::logEmissionProbability(int state, string residue) {
	return logEmissionProbabilities.at(state).at(getEmissionResidueIndex(residue));
}

// double logInitiationProbability(int state)
//  Purpose: 
//		Returns the log of the initiation probability for the state
long double HMMProbabilities::logInitiationProbability(int state) {
	return logInitiationProbabilities[state];
}
	
// double logTransitionProbability(int beginState, int endState)
//  Purpose: 
//		Returns the log of the transition probability for transition from
//		beginState to endState
long double HMMProbabilities::logTransitionProbability(int beginState, int endState) {
	return logTransitionProbabilities[beginState][endState];
}

// long double dSegmentScore(int readStarts)
//  Purpose: 
//		Returns the D-Segment score for the readStarts
long double HMMProbabilities::dSegmentScore(int readStarts) {

	// Get Score contribution form state1
	long double state1Score =
		log(
			 emissionProbabilities.at(1).at(readStarts)
			 * transitionProbability(1, 1)
		) 
		/ log(2);

	// Get Score contribution form state2
	long double state2Score =
		log(
			 emissionProbabilities.at(2).at(readStarts)
			 * transitionProbability(2, 2)
		) 
		/ log(2);

	return state2Score - state1Score; 
}

// setEmissionProbability(int state, char residue, double value)
//  Purpose: 
//		Sets the emission probability for the state and residue to value
//	Postconditions:
//		emissionProbabilites - value set for state/residue
//		logEmissionProbabilites - value set for state/residue
void HMMProbabilities::setEmissionProbability(int state, string residue, long double value) {
	emissionProbabilities[state][getEmissionResidueIndex(residue)] = value;
	double logVal;
	if (value == 0)
		logVal = std::numeric_limits<double>::quiet_NaN();
	else
		logVal = log(value);
	logEmissionProbabilities[state][getEmissionResidueIndex(residue)] = logVal;
}

// setInitiationProbability(int state, double value)
//  Purpose: 
//		Sets the initiation probability for the state to value
//	Postconditions:
//		initiationProbabilites - value set for state
//		logInitiationProbabilites - value set for state
void HMMProbabilities::setInitiationProbability(int state, long double value) {
	initiationProbabilities[state] = value;
	double logVal;
	if (value == 0)
		logVal = std::numeric_limits<double>::quiet_NaN();
	else
		logVal = log(value);
	logInitiationProbabilities[state] = logVal;
}

// setTransitionProbability(int beginState, int endState, double value)
//  Purpose: 
//		Sets the transition probability from the beginState to the endState
//		to value
//	Postconditions:
//		transitionProbabilites - value set for beginState to endState
//		logTransitionProbabilites - value set for beginState to endState
void HMMProbabilities::setTransitionProbability(int beginState, int endState, long double value) {
	transitionProbabilities[beginState][endState] = value;
	double logVal;
	if (value == 0)
		logVal = std::numeric_limits<double>::quiet_NaN();
	else
		logVal = log(value);
	logTransitionProbabilities[beginState][endState] = logVal;
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
string HMMProbabilities::probabilitiesResultsString() {
	stringstream ss;

	// Begin Model
	ss << "      <model type=\"hmm\">\n";

	// States
	ss << statesResultsString();

	// Probabiltiies
	ss << intitiationProbabiltiesResultsString();
	for (int i = 1; i < numStates; i++)
		ss << transitionProbablitiesResultsString(i);
	for (int i = 1; i < numStates; i++)
		ss << emissionProbablitiesResultsString(i);
	
	// End Model
	ss << "      </model>\n";

	return ss.str();
}

// string statesResultsString()
//  Purpose:
//		Returns a string representing the states
//
//		format:
//			<result type="states">
//				<<state1>>,<<state2>>,...
//			</result>
string HMMProbabilities::statesResultsString() {
	stringstream ss;

	// Header 
	ss << "        <states>";

	// States
	for (int i = 1; i < numStates; i++) {
		ss << i ;

		if ( i < numStates - 1)
		   ss << ",";
	}

	// Footer
	ss << "</states>\n";

	return ss.str();
}

// string intitiationProbabiltiesResultsString()
//  Purpose:
//		Returns a string representing the initiation probablities
//
//		format:
//			<result type="initiation_probabilites">
//				<<state>>=<<initiation probability>>,
//			</result>
string HMMProbabilities::intitiationProbabiltiesResultsString() {
	stringstream ss;
	ss.precision(5);

	// Header 
	ss << "        <initial_state_probabilities>";

	// States
	for (int i = 1; i < numStates; i++) {
		ss
			<< i 
			<< "="
			<< initiationProbability(i);

		if ( i < numStates - 1)
		   ss << ",";
	}

	// Footer
	ss << "</initial_state_probabilities>\n";

	return ss.str();
}

// transitionProbablitiesResultsString(int state)
//  Purpose:
//		Returns a string representing the transition probablities for a state
//
//		format:
//			<result type="transition_probabilites" state="<<state>>">
//				<<to state>>=<<transition probability>>,
//			</result>
string HMMProbabilities::transitionProbablitiesResultsString(int state) {
	stringstream ss;
	ss.precision(6);
	ss.scientific;

	// Header 
	ss << "        <transition_probabilities state=\"" << state << "\">";

	// States
	for (int i = 1; i < numStates; i++) {
		ss
			<< i
			<< "="
			<< transitionProbability(state, i);

		if ( i < numStates - 1)
		   ss << ",";
	}

	// Footer
	ss << "</transition_probabilities>\n";

	return ss.str();
}

// emissionProbablitiesResultsString(int state)
//  Purpose:
//		Returns a string representing the emission probablities for a state
//
//		format:
//			<result type="emission_probabilites" state="<<state>>">
//				<<residue>>=<<emission probability>>,
//			</result>
string HMMProbabilities::emissionProbablitiesResultsString(int state) {
	stringstream ss;
	ss.precision(6);
	ss.scientific;

	// Header 
	ss << "        <emission_probabilities state=\"" << state << "\">";

	// Residues
	for (pair<string, int> mapPair : emissionResidueMap) {
		string& residue = mapPair.first;
		ss << residue << "=" << emissionProbability(state, residue) << ",";
	}

	// Footer
	ss << "</emission_probabilities>\n";

	return ss.str();
}

// map<string, int> createEmissionMap()
//  Purpose: 
//		Creates a map of the index location for a nucleotide emission
//		in the emission probabilities array
void HMMProbabilities::createEmissionResidueMap() {

	emissionResidueMap["0"]= 0;
	emissionResidueMap["1"]= 1;
	emissionResidueMap["2"]= 2;
	emissionResidueMap["3"]= 3;
}

// int getIndex(char residue)
//  Purpose: 
//	  Returns the index in the emission probabilities for the residue
int HMMProbabilities::getEmissionResidueIndex(string residue) {
	return emissionResidueMap.at(residue);
}

void HMMProbabilities::populateEmissionProbabilities(int state, double poissonMean) {

	// Set emission for 0 read starts
	double zeroProb = calculatePoissonProbability(poissonMean, 0);
	setEmissionProbability(state, "0", zeroProb);

	// Set emission for 1 read starts
	double oneProb = calculatePoissonProbability(poissonMean, 1);
	setEmissionProbability(state, "1", oneProb);

	// Set emission for 2 read starts
	double twoProb = calculatePoissonProbability(poissonMean, 2);
	setEmissionProbability(state, "2", twoProb);

	// Set emission for 3 or greater read starts
	double threeProb = 1 - (zeroProb + oneProb + twoProb);
	setEmissionProbability(state, "3", threeProb);
}

double HMMProbabilities::calculatePoissonProbability(double mean, int observedValue) {

	return (pow(mean , observedValue) * exp(-mean)) / factorial(observedValue);

}

int HMMProbabilities::factorial(int value) {
	if (value == 0)
		return 1;

	return value * factorial(value - 1);
}


