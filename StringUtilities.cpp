/*
 * StringUtilties.cpp
 *
 *  The StringUtilities object is a container for various operations
 *  that are typically performed on strings but have no std corresponding
 *  std method.
 *
 *  Created on: Jan 21, 2013
 *      Author: tomkolar
 */

#include "StringUtilities.h"
#include <string>
#include <sstream>
#include <vector>
using namespace std;

// Constuctors
// ==============================================
StringUtilities::StringUtilities() {
}

// Destructor
// =============================================
StringUtilities::~StringUtilities() {
}

// Public Class Methods
// =============================================

// vector<string>& split(const string& s, char delim, vector<string>& elems)
//  Purpose: 
//		Split a string into tokens seperated by delim.  The passed in elems
//		vector is populated with the tokens.
//  Postconditions:
//		elems array will be populated with tokens from string
vector<string>& StringUtilities::split(const string& s, char delim, vector<string>& elems) {
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

// string xmlResult(const string& type, const string& value)
//  Purpose: 
//		Returns an XML Result string in the following format:
//			<result type=" <<type>> "> <<value>> <\result>
string StringUtilities::xmlResult(const string& type, const string& value) {
	stringstream ss;
	ss
		<< "    <result type=\"" << type << "\">" << value << "</result>\n";

	return ss.str();
}

// string xmlResult(const string& type, const double value, const int precision)
//  Purpose: 
//		Returns an XML Result string  with the <<value>> set to precision <<precision>>.
//		The return string has the following format:
//			<result type=" <<type>> "> <<value>> <\result>
string StringUtilities::xmlResult(const string& type, const double value, const int precision){
	stringstream ss;
	ss.precision(precision);
	ss
		<< "    <result type=\"" << type << "\">" << value << "</result>\n";

	return ss.str();
}


// string xmlResultFormatted(const string& type, const string& value)
//  Purpose: 
//		Returns an XML Result string in the following format:
//			<result type=" <<type>> ">
//				<<value>>
//			<\result>
string StringUtilities::xmlResultFormatted(const string& type, const string& value) {
	stringstream ss;
	ss
		<< "    <result type=\"" << type << "\">" << "\n"
		<< "      " << value << "\n" 
		<<"    </result>\n";

	return ss.str();
}
