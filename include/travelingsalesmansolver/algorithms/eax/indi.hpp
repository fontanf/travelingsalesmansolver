/*
 * indi.h
 *   created on: April 24, 2013
 * last updated: May 10, 2020
 *       author: Shujia Liu
 */

#ifndef __INDI__
#define __INDI__

#include "travelingsalesmansolver/distances/commons.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


namespace travelingsalesmansolver
{
namespace eax_ga
{

class TIndi{
public:
	TIndi();
	~TIndi();
	void define( int N );
	TIndi& operator = ( const TIndi& src );
	bool operator == (  const TIndi& indi2 ); // checks if two roads are equivalent

	int fN; // the number of cities
	int** fLink; // fLink[i][] is the two adjacent cities of city i
	Distance fEvaluationValue; // the road length of TSP
};


}
}
#endif
