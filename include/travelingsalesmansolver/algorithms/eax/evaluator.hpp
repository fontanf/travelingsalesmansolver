/*
 * evaluator.h
 *   created on: April 24, 2013
 * last updated: May 10, 2020
 *       author: Shujia Liu
 */

#ifndef __EVALUATOR__
#define __EVALUATOR__

#ifndef __INDI__
#include "travelingsalesmansolver/algorithms/eax/indi.hpp"
#endif

#include "travelingsalesmansolver/distances/commons.hpp"

#include <string.h>
#include <assert.h>
#include <vector>
#include <string>
#include <limits>

namespace travelingsalesmansolver
{
namespace eax_ga
{

using namespace std;

/**
 * Distances are not stored/copied here: 'distance(i, j)' is looked up
 * directly on the original (external) 'Distances' object, dispatched once
 * per 'eax()' call via the 'Distances' template parameter, the same way
 * every other algorithm in this library consumes distances.
 */
template <typename Distances>
class TEvaluator{
public:
	TEvaluator( const Distances& distances, int n );
	void doIt( TIndi& indi ); // sets indi.fEvaluationValue
	void writeTo( FILE* fp, TIndi& indi ); // prints out TSP solution
	vector<int> getTour( TIndi& indi ); // returns the 0-indexed tour without going through a file
	bool checkValid(vector<int>& array, Distance value ); // checks if TSP solution is valid

	/** Distance between two cities. */
	inline Distance distance( int vertex_id_1, int vertex_id_2 ) const { return distances_.distance( vertex_id_1, vertex_id_2 ); }

	int fNearNumMax; // the maximum value of the number of nearby points
	vector<vector<int>> fNearCity; // NearCity[i][k] is the k points that with a shortest distance from point i
	int Ncity; // the number of cities
	vector<int> Array; // the index of best solution

private:
	void computeNearCity(); // computes fNearCity from 'distance()'

	const Distances& distances_;
};

template <typename Distances>
TEvaluator<Distances>::TEvaluator( const Distances& distances, int n ):
	fNearNumMax( 50 ),
	Ncity( n ),
	distances_( distances )
{
	fNearCity.clear();
	for (int i = 0; i < Ncity; i++) {
		vector<int> row(fNearNumMax + 1);
		fNearCity.push_back(row);
	}
	computeNearCity();
}

template <typename Distances>
void TEvaluator<Distances>::computeNearCity() {
	vector<int> checkedN(Ncity);
	int ci, j1, j2, j3;
	int cityNum = 0;
	Distance minDis;
	for( ci = 0; ci < Ncity; ++ci ){
		for( j3 = 0; j3 < Ncity; ++j3 ) checkedN[ j3 ] = 0;
		checkedN[ ci ] = 1;
		fNearCity[ ci ][ 0 ] = ci;
		for( j1 = 1; j1 <= fNearNumMax; ++j1 ) {
			minDis = std::numeric_limits<Distance>::max();
			for( j2 = 0; j2 < Ncity; ++j2 ){
				if( distance( ci, j2 ) <= minDis && checkedN[ j2 ] == 0 ){
					cityNum = j2;
					minDis = distance( ci, j2 );
				}
			}
			fNearCity[ ci ][ j1 ] = cityNum;
			checkedN[ cityNum ] = 1;
		}
	}
}

template <typename Distances>
void TEvaluator<Distances>::doIt( TIndi& indi ) {
	Distance d = 0;
	for( int i = 0; i < Ncity; ++i ) d += distance( i, indi.fLink[i][0] ) + distance( i, indi.fLink[i][1] );
	indi.fEvaluationValue = d/2;
}

template <typename Distances>
vector<int> TEvaluator<Distances>::getTour( TIndi& indi ) {
	vector<int> tour(Ncity);
	int curr = 0, st = 0, count = 0, pre = -1, next;
	while( 1 ) {
		tour[ count++ ] = curr;
		if( count > Ncity ) {
			printf( "Invalid\n" );
			break;
		}
		if( indi.fLink[ curr ][ 0 ] == pre ) next = indi.fLink[ curr ][ 1 ];
		else next = indi.fLink[ curr ][ 0 ];

		pre = curr;
		curr = next;
		if( curr == st ) break;
	}
	return tour;
}

template <typename Distances>
void TEvaluator<Distances>::writeTo( FILE* fp, TIndi& indi ){
	Array.resize(Ncity);
	int curr=0, st=0, count=0, pre=-1, next;
	while( 1 ){
		Array[ count++ ] = curr + 1;
		if( count > Ncity ){
			printf( "Invalid\n" );
			return;
		}
		if( indi.fLink[ curr ][ 0 ] == pre ) next = indi.fLink[ curr ][ 1 ];
		else next = indi.fLink[ curr ][ 0 ];

		pre = curr;
		curr = next;
		if( curr == st ) break;
	}
	if( this->checkValid( Array, indi.fEvaluationValue ) == false )
		printf( "Individual is invalid \n" );

	fprintf( fp, "%d %lld\n", indi.fN, (long long)indi.fEvaluationValue );
	for( int i = 0; i < indi.fN; ++i )
		fprintf( fp, "%d ", Array[ i ] );
	fprintf( fp, "\n" );
}

template <typename Distances>
bool TEvaluator<Distances>::checkValid(vector<int>& array, Distance value) {
	int *check=new int[Ncity];
	for( int i = 0; i < Ncity; ++i ) check[ i ] = 0;
	for( int i = 0; i < Ncity; ++i ) ++check[ array[ i ]-1 ];
	for( int i = 0; i < Ncity; ++i )
		if( check[ i ] != 1 ) return false;
	Distance total_distance = 0;
	for( int i = 0; i < Ncity-1; ++i )
		total_distance += distance( array[ i ]-1, array[ i+1 ]-1 );

	total_distance += distance( array[ Ncity-1 ]-1, array[ 0 ]-1 );

	delete [] check;
	if( total_distance != value ) return false;
	return true;
}

}
}
#endif
