
#ifndef VEC3F_H
#define VEC3F_H

#include <iostream>

// DUMMY Vec3f class.

struct Vec3f
{
	Vec3f( float x = 0.f, float y = 0.f, float z = 0.f ) :
		m_x( x ),
		m_y( y ),
		m_z( z )
	{}

	float m_x, m_y, m_z;

	friend std::istream & operator >> ( std::istream & in, Vec3f & vec )
	{
		float x, y, z;

		if( in >> x >> y >> z )
		{
			vec.m_x = x;
			vec.m_y = y;
			vec.m_z = z;
		}
		
		return in;
	}

	friend std::ostream & operator << ( std::ostream & out, Vec3f const & vec )
	{
		return out << vec.m_x << " " << vec.m_y << " " << vec.m_z;
	}
};

#endif
