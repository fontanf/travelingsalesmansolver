/*
 * environment.h
 *   created on: April 24, 2013
 * last updated: May 10, 2020
 *       author: Shujia Liu
 */

#ifndef __ENVIRONMENT__
#define __ENVIRONMENT__

#ifndef __INDI__
#include "travelingsalesmansolver/algorithms/eax/indi.hpp"
#endif

#ifndef __RAND__
#include "travelingsalesmansolver/algorithms/eax/randomize.hpp"
#endif

#ifndef __EVALUATOR__
#include "travelingsalesmansolver/algorithms/eax/evaluator.hpp"
#endif

#ifndef __Cross__
#include "travelingsalesmansolver/algorithms/eax/cross.hpp"
#endif

#ifndef __KOPT__
#include "travelingsalesmansolver/algorithms/eax/kopt.hpp"
#endif

#include "travelingsalesmansolver/distances/commons.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <vector>


namespace travelingsalesmansolver
{
namespace eax_ga
{

/**
 * 'distances' is looked up directly on the original (external) 'Distances'
 * object throughout (see 'TEvaluator'); nothing here ever copies it into a
 * separate matrix.
 */
template <typename Distances>
class TEnvironment{
public:
	TEnvironment( const Distances& distances, int n );
	~TEnvironment();

	void define(); // global initialization
	void doIt(); // entry point of genetic algorithm
	void init(); // initializes genetic algorithm
	bool terminationCondition(); // condition to termination the genetic algorithm
	void setAverageBest(); // sets the average distance of TSP and the shortest distance of TSP

	void initPop(); // initializes population
	void selectForMating(); // selects parents
	void generateKids( int s ); // generates and selects children
	void getEdgeFreq(); // calculates the frequency of each edge

	void printOn( int n ); // logs out results
	void writeBest(); // logs out the best TSP solution

	TEvaluator<Distances>* fEvaluator; // distance of each edge
	TCross<Distances>* tCross; // intersection of edge sets
	TKopt<Distances>* tKopt; // local search (2-opt neighborhood)

	int Npop; // the number of population
	int Nch; // the number of children generated from one pair of parent
	TIndi* tCurPop; // the member of the current population
	TIndi tBest; // the best solution from the current population
	int fCurNumOfGen; // the number of generations of the current population
	long int fAccumurateNumCh; // accumulated number of children so far

	int fBestNumOfGen; // the number of generations that generates the current best solution
	long int fBestAccumeratedNumCh; // accumulated number of generations of the current best solution
	vector<vector<int>> fEdgeFreq; // edge frequency of a population
	double fAverageValue; // average road length of TSP in a population
	Distance fBestValue; // road length of the best solution in a population
	int fBestIndex;	// index of the best solution in a population

	vector<int> fIndexForMating; // list for edge cross operation
	int fStagBest; // accumulated number of generations that doesn't generate a better solution compared to previous generation
	int fFlagC[ 10 ]; // EAX method and selection strategy
	int fStage; // the current step of genetic algorithm
	int fMaxStagBest; // the genetic algorithm goes into the next stage when fStagBest == fMaxStagBest
	int fCurNumOfGen1; // the number of generations when stage 1 completes

	clock_t fTimeStart, fTimeInit, fTimeEnd; // saves calculation time
	double fTimeLimit; // wall-clock time limit in seconds, negative means no limit
};

template <typename Distances>
TEnvironment<Distances>::TEnvironment( const Distances& distances, int n ){
	fEvaluator = new TEvaluator<Distances>( distances, n );
	fTimeLimit = -1.0;
}

template <typename Distances>
TEnvironment<Distances>::~TEnvironment(){
	delete [] tCurPop;
	delete fEvaluator;
	delete tCross;
	delete tKopt;
}

template <typename Distances>
void TEnvironment<Distances>::define(){
	int N = fEvaluator->Ncity;
	fIndexForMating.resize(Npop + 1);

	tCurPop = new TIndi [ Npop ];
	for ( int i = 0; i < Npop; ++i ) {
		tCurPop[i].define( N );
	}
	tBest.define( N );
	tCross = new TCross<Distances>( N );
	tCross->eval = fEvaluator;
	tCross->Npop = Npop;
	tKopt = new TKopt<Distances>( N );
	tKopt->eval = fEvaluator;
	tKopt->setInvNearList();

	fEdgeFreq.clear();
	for (int i = 0; i < N; i++) {
		vector<int> row(N);
		fEdgeFreq.push_back(row);
	}
}

template <typename Distances>
void TEnvironment<Distances>::doIt(){
	this->fTimeStart = clock();
	this->initPop(); // initializes population
	this->fTimeInit = clock();
	this->init();

	this->getEdgeFreq();
	while( 1 ){
		this->setAverageBest();
		if( this->terminationCondition() ) break;

		this->selectForMating();
		for( int s =0; s < Npop; ++s ) this->generateKids( s );

		++fCurNumOfGen;
	}
	this->fTimeEnd = clock();
}

template <typename Distances>
void TEnvironment<Distances>::init(){
	fAccumurateNumCh = 0;
	fCurNumOfGen = 0;
	fStagBest = 0;
	fMaxStagBest = 0;
	fStage = 1; // sets stage to 1
	fFlagC[ 0 ] = 4; // maintains population diversity	1:Greedy, 2:---, 3:Distance, 4:Entropy
	fFlagC[ 1 ] = 1; // the type of Eset: 1:Single-AB, 2:Block2
}

template <typename Distances>
bool TEnvironment<Distances>::terminationCondition(){
	if ( fTimeLimit >= 0.0
			&& (double)(clock() - fTimeStart) / CLOCKS_PER_SEC >= fTimeLimit )
		return true;
	if ( fAverageValue - fBestValue < 0.001 )  return true;
	if( fStage == 1 ){
		if( fStagBest == int(1500/Nch) && fMaxStagBest == 0 ) // 1500/Nch
			fMaxStagBest =int( fCurNumOfGen / 10 ); // fMaxStagBest = G/10
		else if( fMaxStagBest != 0 && fMaxStagBest <= fStagBest ){
			fStagBest = 0;
			fMaxStagBest = 0;
			fCurNumOfGen1 = fCurNumOfGen;
			fFlagC[ 1 ] = 2;
			fStage = 2;
		}
		return false;
	}
	if( fStage == 2 ){
		if( fStagBest == int(1500/Nch) && fMaxStagBest == 0 )			// 1500/Nch
			fMaxStagBest = int( (fCurNumOfGen - fCurNumOfGen1) / 10 );	// fMaxStagBest = G/10
		else if( fMaxStagBest != 0 && fMaxStagBest <= fStagBest ) return true;
		return false;
	}

	return true;
}

template <typename Distances>
void TEnvironment<Distances>::setAverageBest(){
	Distance stockBest = tBest.fEvaluationValue;
	fAverageValue = 0.0;
	fBestIndex = 0;
	fBestValue = tCurPop[0].fEvaluationValue;
	for(int i = 0; i < Npop; ++i ){
		fAverageValue += tCurPop[i].fEvaluationValue;
		if( tCurPop[i].fEvaluationValue < fBestValue ){
			fBestIndex = i;
			fBestValue = tCurPop[i].fEvaluationValue;
		}
	}
	tBest = tCurPop[ fBestIndex ];
	fAverageValue /= (double)Npop;
	if( tBest.fEvaluationValue < stockBest ){
		fStagBest = 0;
		fBestNumOfGen = fCurNumOfGen;
		fBestAccumeratedNumCh = fAccumurateNumCh;
	}
	else ++fStagBest;
}

template <typename Distances>
void TEnvironment<Distances>::initPop(){
	for ( int i = 0; i < Npop; ++i ){
		tKopt->makeRandSol( tCurPop[ i ] ); // randomly sets a route
		tKopt->doIt( tCurPop[ i ] ); // local search (2-opt neighborhood)
	}
}

template <typename Distances>
void TEnvironment<Distances>::selectForMating(){
	tRand->permutation( fIndexForMating, Npop, Npop );
	fIndexForMating[ Npop ] = fIndexForMating[ 0 ];
}

template <typename Distances>
void TEnvironment<Distances>::generateKids( int s ){
	// tCurPop[fIndexForMating[s]] gets replaced by the best solution in tCross->DoIt()
	// fEdgeFreq[][] gets updated at the same time
	tCross->setParents( tCurPop[fIndexForMating[s]], tCurPop[fIndexForMating[s+1]], fFlagC, Nch );
	tCross->doIt( tCurPop[fIndexForMating[s]], tCurPop[fIndexForMating[s+1]], Nch, 1, fFlagC, fEdgeFreq );
	fAccumurateNumCh += tCross->fNumOfGeneratedCh;
}

template <typename Distances>
void TEnvironment<Distances>::getEdgeFreq(){
	int  k0, k1, N = fEvaluator->Ncity;
	for( int j1 = 0; j1 < N; ++j1 )
		for( int j2 = 0; j2 < N; ++j2 )
			fEdgeFreq[ j1 ][ j2 ] = 0;

	for( int i = 0; i < Npop; ++i )
		for(int j = 0; j < N; ++j ){
			k0 = tCurPop[ i ].fLink[ j ][ 0 ];
			k1 = tCurPop[ i ].fLink[ j ][ 1 ];
			++fEdgeFreq[ j ][ k0 ];
			++fEdgeFreq[ j ][ k1 ];
		}
}

template <typename Distances>
void TEnvironment<Distances>::printOn( int n ){
	printf( "n = %d val = %lld Gen = %d Time = %d %d\n" , n, (long long)tBest.fEvaluationValue, fCurNumOfGen,
		(int)((double)(this->fTimeInit - this->fTimeStart)/(double)CLOCKS_PER_SEC),
		(int)((double)(this->fTimeEnd - this->fTimeStart)/(double)CLOCKS_PER_SEC) );
	fflush(stdout);

}

template <typename Distances>
void TEnvironment<Distances>::writeBest(){
	FILE *fp;
	char filename[ 80 ];

	sprintf( filename, "bestSolution.txt" );
	fp = fopen( filename, "a");
	fEvaluator->writeTo( fp, tBest );
	fclose( fp );
}

}
}
#endif
