
#ifndef VEC3F_FILEIO_H
#define VEC3F_FILEIO_H

#include <iostream>
#include <string>
#include <cstddef>
#include <vector>
#include <sstream>
#include <iterator>
#include <fstream>
#include <stdexcept>

#include "Vec3f.h"

typedef std::vector< Vec3f > VectorContainerVec3f;

// i) Reads in file (e.g. text file (.txt) or contour file (.con)).
// ii) Parses it line by line.
// iii) Extracts the first three floating point values (ignoring all other values)
// iv) Stores these values as a Vec3f in a vector container
void loadVec3fFromFile( VectorContainerVec3f & vecs, std::string const & fileName )
{
	using std::string;
	using std::stringstream;
	using std::istream_iterator;

	std::ifstream file( fileName );

	if( !file )
	{
		throw std::runtime_error( "Unable to open file." );
	}

	string line;
	size_t index;
	stringstream ss( std::ios_base::in );

	size_t lineNum = 0;
	vecs.clear();

	while( getline( file, line ) )
	{	
		++lineNum;
		
		// remove comments	
		index = line.find_first_of( "#" );
		if( index != string::npos )
		{
			line.erase( index, string::npos );
		}

		// removes leading/tailing junk
		line.erase( 0, line.find_first_not_of( " \t\r\n\v\f" ) );
		index = line.find_last_not_of( " \t\r\n\v\f" )+1;
		if( index != string::npos )
		{
			line.erase( index, string::npos );
		}
		
		if( line.empty() )
		{
			continue; // empty or commented out line
		}
		
		ss.str( line );
		ss.clear();

		try{
			Vec3f v;
			ss >> v;

			if( !ss )
			{
			/*	throw std::runtime_error( "Error read file: "
							 + line
							 + " (line: "
							 + std::to_string(lineNum)
							 +  ")" );
			*/}	
			else
			{
				vecs.push_back( v );
			}
		} catch (const std::exception & e) { std::cerr << e.what() << std::endl; }
	}
	file.close();
}

#endif
