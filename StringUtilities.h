/*
 * StringUtilties.h
 *
 *  The StringUtilities object is a container for various operations
 *  that are typically performed on strings but have no std corresponding
 *  std method.
 *
 *  Created on: Jan 21, 2013
 *      Author: tomkolar
 */

#ifndef STRINGUTILITIES_H_
#define STRINGUTILITIES_H_

#include <string>
#include <vector>
using namespace std;

class StringUtilities {
public:

	// Constuctors
	// ==============================================
	StringUtilities();

	// Destructor
	// =============================================
	virtual ~StringUtilities();

	// Public Class Methods
	// =============================================

	// vector<string>& split(const string& s, char delim, vector<string>& elems)
	//  Purpose: 
	//		Split a string into tokens seperated by delim.  The passed in elems
	//		vector is populated with the tokens.
	//  Postconditions:
	//		elems array will be populated with tokens from string.
	static vector<string>& split(const string& s, char delim, vector<string>& elems);

	// string xmlResult(const string& type, const string& value)
	//  Purpose: 
	//		Returns an XML Result string in the following format:
	//			<result type=" <<type>> "> <<value>> <\result>
	static string xmlResult(const string& type, const string& value);

	// string xmlResult(const string& type, const double value, const int precision)
	//  Purpose: 
	//		Returns an XML Result string  with the <<value>> set to precision <<precision>>.
	//		The return string has the following format:
	//			<result type=" <<type>> "> <<value>> <\result>
	static string xmlResult(const string& type, const double value, const int precision);

	// string xmlResultFormatted(const string& type, const string& value)
	//  Purpose: 
	//		Returns an XML Result string in the following format:
	//			<result type=" <<type>> ">
	//				<<value>>
	//			<\result>
	static string xmlResultFormatted(const string& type, const string& value);
};

#endif /* STRINGUTILITIES_H_ */
