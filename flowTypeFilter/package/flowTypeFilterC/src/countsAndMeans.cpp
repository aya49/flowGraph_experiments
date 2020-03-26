// [[Rcpp::depends(BH)]]

#include "numberGenerator.h"
#include "exceptions.h"
#include <vector>
#include <map>
#include <boost/dynamic_bitset.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include "countsAndMeans.h" 
#include <R_ext/Print.h>
#include <Rmath.h>

using namespace std;
using namespace tfl;
using namespace boost;
using namespace Rcpp ;
namespace bg = boost::geometry;



/*****************************************************************************
 * Utility code to allow sorting of data for MFI determination while  
 * to cout
 ****************************************************************************/
typedef std::pair<double,int> valueWithIndex;

bool comparator ( const valueWithIndex& l, const valueWithIndex& r)
   { return l.first < r.first; }



/*****************************************************************************
 * Determine whether an event falls within a given range 
 *
 ****************************************************************************/

bool is_in_range(const double data, const vector<double> &thresholds, const int segment, const int segment_count)
{
	if (segment == 1)
	{
		
		return data <= thresholds[0];
	}
	else if (segment == segment_count)
	{
		
		return data > thresholds[segment_count - 2];
	}
	else
	{
		return data <= thresholds[segment - 1] && data > thresholds[segment - 2];
	}
	throw "can not happen!";
}

/*****************************************************************************
 * Determine whether an event falls within a given filter
 *
 ****************************************************************************/

bool is_in_filter(const double data_x, const double data_y, const vector<double> poly_points , const int segment)
{

    typedef bg::model::d2::point_xy<double> point_type;
    typedef bg::model::polygon<point_type> polygon_type;

    point_type p(data_x, data_y);

    //convert the vector of double from threshold into a vector of points
    std::vector<point_type> points;
    for (int i = 0; i < poly_points.size()/2; i++)
    {
      // poly_point contains 2 columns of data
      // point_type(x, y)
      points.push_back(point_type(poly_points[i], poly_points[poly_points.size()/2+i]));
    }

    // use vector of points to create the filter polygon
    polygon_type poly;
    bg::assign_points(poly, points);
    
    //use function "covered_by" from boost library for segmenting
    bool in_poly = bg::covered_by(p, poly);
    
    if (segment == 1)
    {
        return !in_poly;
    }
    else if (segment == 2)
    {
        return in_poly;
    }
    else if (segment >= 2)
    {
        throw "Cannot split data to more than 2, when polygone filter is provided.";
    }
    throw "can not happen!";
}


/*****************************************************************************
 * Generating bitstring specifying which events in phenotype
 *
 ****************************************************************************/



dynamic_bitset<> terminal_phenotype(
    const char *phenotype, //
    const double *data, //List mode flow cytometry data
    vector<vector<double> > &thresholds, //Thresholds, variable length vector for each marker 
    const int gate_ind, // First gate index in threshold list (value -1 if there's no gates)
    const int *thresChannels, //lables of marker channels for each threshold
    const int data_count, //Number of events (rows in list mode data)
    const int marker_count, //Number of markers (columns in list mode data)
    const int *base // Array of number of thresholds per marker (length = marker_count)
    )
{
	dynamic_bitset<> result(data_count);
	int marker = -1;
  for (register int i = 0; i < marker_count; i++)
	{
		if (phenotype[i])
		{
			marker = i;
			break;
		}
	}

	if (marker == -1)
		throw "no marker present!!!";

  int marker_value = phenotype[marker];
  
  if (gate_ind == -1) //only thresholds
  {
    //determine which threshold segment are the events in (smaller or bigger than threshold)
    for (register int i = 0; i < data_count; i++)
    {
      // marker columns in data (extrated by PropMarkers) are in the same order as in thresholds
      result[i] = (is_in_range(data[marker * data_count + i], thresholds[marker], marker_value, base[marker] - 1) ? 1 : 0);
    }
  }
  else //thresholds and filters
  {	
    if (marker < gate_ind) // this marker use threshold
  	{
  		for (register int j = 0; j < data_count; j++)
  		{
  		  int col = thresChannels[marker]-1; //which column of data should we extract
  		  result[j] = (is_in_range(data[col * data_count + j], thresholds[marker], marker_value, base[marker] - 1) ? 1 : 0);
  		}
  	}
  	else // this marker use filter
  	{
  		//determine which gate segment are the events in (in or out of the polygon)
  		int gate_id = marker - gate_ind;  //which gate it is
  		for (register int j = 0; j < data_count; j++)
  		{
  		  //which columns of data should we extract
  		  int col_x = thresChannels[gate_ind + gate_id * 2] -1; 
  		  int col_y = thresChannels[gate_ind + gate_id * 2 + 1] -1; 
  		  
  		  result[j] = (is_in_filter(data[col_x * data_count + j], //gat point.x in data according to threChannels
  		            data[col_y * data_count + j], //get point.y in data according to threChannels
  		            thresholds[marker], //vector in this gate arrays as points of polygon (poly_x and poly_y)
  		            marker_value) ? 1 : 0);
  	  }
  	 }
  	}

	return result;
}

/*****************************************************************************
 * 
 *
 ****************************************************************************/


void elements(
    //Inputs:
    NumberGenerator &n, //Number generator for sequence of phenotypes
    //Outputs:
    char *parent, //
    char *terminal //
    )
{
	strcpy(parent, n.neighbor(0).text());
	for (int i = 0; i < n.getLength(); i++)
	{
		terminal[i] = '0';
	}
	int dif_idx = n.firstDiff(n.neighbor(0));
	terminal[dif_idx] = n.text()[dif_idx];
	terminal[n.getLength()] = 0;
}


/*****************************************************************************
 * Calculate Median Fluorescence Intensities for the given phenotype 
 * 
 ****************************************************************************/

void calculateMFIs(
        const dynamic_bitset<> &phenotype, //bitstring specifying which events in phenotype
        const vector< vector<double>* > &sorted_channels, //Channel data, not actually sorted in order
        const vector< map<int,int>* > &event_map, //Map from sorted order to bitstring order
        NumericVector &return_mfis
	)
{
  // reminder: I construct mfis here, should be destructed outside. -- not any more :)

    //2. For each MFI channel: 
    // 2.a. Create blank vector for event values
    // 2.b. Loop over bitset, pushing values into blank vector when 1s
    // 2.c. If odd no/events, take middle value; if even, 
    // take sum(middle, middle+1)/2   [NOte: use count on bitset to find length]

    int count = phenotype.count();
    // Only compute MFI if the phenotype is non-empty: 
    if(count)
    {
        int midpoint = count / 2;
    

        for(size_t channel=0; channel < sorted_channels.size(); channel++)
        {
            vector<double> sorted_subset;
            vector<double> *this_channel = sorted_channels[channel];
            map<int, int> *this_map = event_map[channel];
            
	    //Take subset of sorted channel based on phenotype's set bits:
            for (size_t event=0; event < phenotype.size(); event++)
            {
	      int sorted_event = (*this_map)[event];

                if(phenotype[sorted_event])
		{
		  sorted_subset.push_back((*this_channel)[sorted_event]);
		}
            }
	    //Either take the exact midpoint if odd length, or average of mid two points if even:
            if (count % 2)
            	return_mfis[channel] = sorted_subset[midpoint];
                
            else
            	return_mfis[channel] = (sorted_subset[midpoint] + sorted_subset[midpoint-1]) /2.0;
                
        }
    }
    //If subset is of size 0, set mfi = NA for all channels
    else
    {
        for(size_t channel=0; channel < sorted_channels.size(); channel++)
        	return_mfis[channel] = NA_REAL;
    }

    //return mfis;

}


/*****************************************************************************
 * Calculate cell counts and MFIs for all populations specified
 *
 * 
 ****************************************************************************/


void calculateCountsAndMeans(
                        const int marker_count, //total number of markers to threshold
                        const int *base, //number of thresholds for each channel 
                        const int data_count, //number of rows in list mode flow data 
                        const int max_level, //only generate phenotypes of up to this many markers
                        vector<vector<double> > &thresholds, //values of thresholds defining partitions, given per channel
                        const int gateInd, //Index of first gate filter in threshold vectors, be -1 if there's no gates
                        const int *thresChannels, //lables of marker channels for each threshold
                        double *data, //"list mode" flow data (matrix of expression values for cells by channel)
                        NumericMatrix mfi_data, //"list mode" flow data, only MFI channels
                        const int verbose,

                        //Outputs:
                        IntegerVector &counts, // 
                        NumericMatrix &mfis, // 
                        StringVector &pheno_codes //byte codes representing the immunophenotypes

                         )

{
        //Locals:
	map<string, dynamic_bitset<> > firstLevel; //To store all phenotypes made up of just one marker
	map<string, dynamic_bitset<> > currentLevel; //To store all phenotypes in the current level (n markers)
	map<string, dynamic_bitset<> > previousLevel; //To store all phenotypes in the previous level (n-1 markers)

	 char p1[1000], p2[1000]; //Space to store the two phenotypes currently being worked on
	 //NOTE: If somehow more than 1000 markers were used, this would overflow!
	 string currPheno;
	 unsigned int popCounter = 0; //Counter keeping track of current phenotype (for indexing outputs)
	 //Note: if more than 2^32 populations, this will cause errors
	 NumericVector current_mfis(mfi_data.ncol()); // Temporary storags for MFIs returned from calculateMFIs
	 

	 //Set up pre-sorted MFI channels for fast computation of MFI:

	 vector< vector<double>* > sorted_channels; //Note: not actually sorted; sorting is stored in index_map
	 vector< map<int,int>* > index_map;
	 if(mfi_data.nrow()==data_count)
	 {
	     for(int i=0; i < mfi_data.ncol(); i++)
	     {
	       vector<double> *channel = new vector<double>(mfi_data(_,i).begin(), mfi_data(_,i).begin() + mfi_data.nrow());
		     sorted_channels.push_back(channel);
		     vector<valueWithIndex> *channel_index = new vector<valueWithIndex>();
		    
		    //Store the channel in a paired vector, with the index for each value: 
		     for (int j=0; j< mfi_data.nrow(); j++)
		     {
		       valueWithIndex this_value((*channel)[j], j);
			(*channel_index).push_back(this_value);
		     }

		    //Sort the channel by value, but keep indices:
		    sort(channel_index->begin(), channel_index->end());

		     //Now find where the indices ended up, and store these in the map:
		    map<int,int> *channel_map = new map<int, int>();
		     for (int j=0; j< mfi_data.nrow(); j++)
		     {
		       int this_index = (*channel_index)[j].second;
		       (*channel_map)[j] = this_index;
		     }

		     index_map.push_back(channel_map);
		     delete channel_index;
	     }
	 }

	//Add root phenotype.
	string root_phenotype = string(marker_count, '0');
	pheno_codes(popCounter) = root_phenotype;
	counts[popCounter] = data_count;
	dynamic_bitset<> all_data_points(data_count);
	all_data_points.set();
	
	calculateMFIs(all_data_points, sorted_channels, index_map, current_mfis);
	mfis(popCounter,_) = current_mfis;
	
	popCounter++;

	 //Generate all first level (terminal / single-marker) phenotypes:

	NumberGenerator n(marker_count, base, 1);
	n.first();
	while(n.hasNext())
	{
		n.next();
        firstLevel[n.text()] = terminal_phenotype(n.data(), data, thresholds, gateInd, thresChannels, data_count, marker_count, base);

		//Add counts, mfis and phenotype codes to return containers:
		counts[popCounter] = firstLevel[n.text()].count();

		calculateMFIs(firstLevel[n.text()], sorted_channels, index_map, current_mfis);
		mfis(popCounter,_) = current_mfis;

		currPheno = string(n.text());
		pheno_codes(popCounter) = currPheno; 
		popCounter++;
                
	}
	//---------------------------------------------------------------------------

	if(verbose)
	  Rprintf("First level done. Identified phenotypes: %i\n", popCounter);

	n = NumberGenerator(marker_count, base, 2);
	n.first();
	while(n.hasNext())
	{
		n.next();
		//cerr << n.text() << " " << n.neighbor(0).text() << " " << n.text()[n.firstDiff(n.neighbor(0))] << endl;
		elements(n, p1, p2);
		//cerr << n.text() << " " << p1 << " " << p2 << endl;
		currentLevel[n.text()] = firstLevel[string(p1)] & firstLevel[string(p2)];
        counts[popCounter] = currentLevel[n.text()].count();
                
        calculateMFIs(currentLevel[n.text()], sorted_channels, index_map, current_mfis);
		mfis(popCounter,_) = current_mfis;
        currPheno = string(n.text());
        pheno_codes(popCounter) = currPheno; 

        popCounter++;
	
	}
        
	//dump(firstLevel);
        //---------------------------------------------------------------------------

	if (verbose)
	  Rprintf("Second level done. Identified phenotypes: %i\n", popCounter);
        //Dynamic programming step:
        //Generate all phenotypes in level i, using previous two levels
	for (int i = 3; i < max_level + 1; i++)
	{
	  //Transfer currentLevel to previousLevel and clear currentLevel for new level to be inserted: 
                previousLevel.insert(currentLevel.begin(), currentLevel.end());
                currentLevel.clear();

		NumberGenerator n(marker_count, base, i);
		n.first();
                //Generate all phenotypes for level i:
		while (n.hasNext())
		{
			n.next();
			elements(n, p1, p2);
			currentLevel[n.text()] = previousLevel[string(p1)] & firstLevel[string(p2)];
                        
                        counts[popCounter] = currentLevel[n.text()].count();
                        
                        calculateMFIs(currentLevel[n.text()], sorted_channels, index_map, current_mfis);
				mfis(popCounter,_) = current_mfis;

                        currPheno = string(n.text());
                        pheno_codes(popCounter) = currPheno; 
                        popCounter++;
		}
       
		if( verbose)
	Rprintf("Level %i done. Identified phenotypes: %i\n", i, popCounter);
	}
        //---------------------------------------------------------------------------
        if(verbose)
	Rprintf("All levels done. Identified phenotypes: %i\n", popCounter);
}


