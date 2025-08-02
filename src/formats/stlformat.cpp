//
// STL Plugin for Open Babel
// Copyright (C) 2015 M J Harvey,
// Acellera Ltd
// m.j.harvey ( at ) acellera.com
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA  02110-1301, USA.
//
// $Author$
// $Date$
// $Revision$
//


#include <openbabel/math/vector3.h>
#include <openbabel/mol.h>
#include <openbabel/atom.h>
#include <openbabel/elements.h>

#include <openbabel/obconversion.h>
#include <openbabel/obmolecformat.h>

#include <openbabel/obiter.h>
#include <openbabel/data.h>

#include <iostream>
#include <ostream>
#include <vector>
#include <array>

#include <cstdlib>
#include <cmath>

// uint8_t and uint16_t are not defined for earlier versions of msvc
#if defined(_MSC_VER) && _MSC_VER <= 1600
  typedef unsigned __int8 uint8_t;
  typedef unsigned __int16 uint16_t;
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


using namespace std;

namespace OpenBabel
{


  //==============================================================================
  /// Class to output a point cloud on a surface around a molecule

  class STLFormat : public OpenBabel::OBMoleculeFormat
  {
    public:
      // Non-standard STL extension for colour: 5 bits per col, bit 15 set
      template<uint8_t Red, uint8_t Green, uint8_t Blue> 
      static constexpr uint16_t colour() {
        return 0x800 | ((Red & 0x1F) << 10) | ((Green & 0x1F) << 5) | (Blue & 0x1F);
      };

    public:
      STLFormat()
      {
        OpenBabel::OBConversion::RegisterFormat( "stl", this );
        /*
        OBConversion::RegisterOptionParam("p", this, 1, OBConversion::OUTOPTIONS);
        OBConversion::RegisterOptionParam("s", this, 1, OBConversion::OUTOPTIONS);
        OBConversion::RegisterOptionParam("c", this, 1, OBConversion::OUTOPTIONS);
        */
      }

      /// Return description.
      const char* Description() override  // required
      {
        return
          "STL 3D-printing format\n"
          "The STereoLithography format developed by 3D Systems\n\n"
          "Write Options, e.g. -xc\n"
          "  p <radius> radius for probe particle (default 0.0 A)\n"
          "  s <scale> scale-factor for VDW radius (default 1.0 A)\n"
          "  c add CPK colours\n\n";
      }

      const char* SpecificationURL() override
      {
        return "http://www.ennex.com/~fabbers/StL.asp";
      }

      /// Return MIME type, NULL in this case.
      const char* GetMIMEType() override { return "application/sla"; }

      /// Return read/write flag: read only.
      unsigned int Flags() override
      {
        return WRITEONEONLY | NOTREADABLE;
      }

      /// Skip to object: used for multi-object file formats.
      int SkipObjects(int n, OpenBabel::OBConversion* pConv) override { return 0; }

      /// Read: always return false.
      bool ReadMolecule(OpenBabel::OBBase*, OpenBabel::OBConversion*) override
      {
        return false;
      }

      /// Write.
      bool WriteMolecule(OpenBabel::OBBase*, OpenBabel::OBConversion*) override;
  };


  // Global variable used to register STL format.
  STLFormat theSTLFormat;


  // Produces triangles over the surface of a sphere


  struct Triangle {
    vector3 a, b, c;
    uint16_t col;
    Triangle( vector3 x, vector3 y, vector3 z, uint16_t colour ) : a{x}, b{y}, c{z}, col{colour} {
    }
  };

  static void map_sphere ( vector<Triangle> &triangles, vector3 origin, double r, uint16_t col )
  {
    static constexpr std::size_t LongitudeSteps{ 144 };
    static constexpr std::size_t LatitudeSteps{LongitudeSteps / 2};

    vector<vector3> points;

    double theta =  ( 2*M_PI / LongitudeSteps );
    int p2 = LongitudeSteps / 2;
    int r2 = LatitudeSteps / 2;
    for(int y = -r2; y < r2; ++y) {
      double cy = cos(y*theta);
      double cy1 = cos((y+1)*theta);
      double sy = sin(y*theta);
      double sy1 = sin((y+1)*theta);

      for(int i = -p2; i < p2; ++i) {
        double ci = cos(i*theta);
        double si = sin(i*theta);
        points.emplace_back(origin[0] + r * ci*cy , origin[1] + r*sy , origin[2] + r * si*cy);
        points.emplace_back(origin[0] + r * ci*cy1, origin[1] + r*sy1, origin[2] + r * si*cy1);
      }
    }

    for( int i=0 ; i < points.size()-2; i++ ) {
      // Order the points to obey the 'right-hand rule', so normals all face outwards
      if( i%2 == 0 ) {
        triangles.emplace_back(points[i], points[i+1], points[i+2], col);
      }
      else {
        triangles.emplace_back(points[i+2], points[i+1], points[i], col);
      }
    }
  }

  static uint16_t cpk_colour( unsigned int atomicnum ) {
    // CPK colouring
    switch( atomicnum ) {
      case 1:  return STLFormat::colour<0x1F, 0x1F, 0x1F>(); // H, WHITE
      case 6:  return STLFormat::colour<0x04, 0x04, 0x04>(); // C, (almost) BLACK
      case 7:  return STLFormat::colour<0x12, 0x1A, 0x1F>(); // N, SKY BLUE
      case 8:  return STLFormat::colour<0x1F, 0x00, 0x00>(); // O, RED
      case 16: return STLFormat::colour<0x1F, 0x1F, 0x00>(); // S, YELLOW
      case 15: return STLFormat::colour<0x1F, 0x00, 0x1F>(); // P, PURPLE
      case 9:  return STLFormat::colour<0x00, 0x1F, 0x00>(); // F, LIGHT GREEN
      case 17: return STLFormat::colour<0x00, 0x17, 0x00>(); // Cl, MEDIUM GREEN
      case 35: return STLFormat::colour<0x00, 0x0F, 0x00>(); // Br, MEDIUM DARK GREEN
      case 53: return STLFormat::colour<0x00, 0x07, 0x00>(); // I, DARK GREEN
      case 26: // Fe
      case 27: // Co
      case 28: // Ni
      case 29:	return STLFormat::colour<0x18, 0x18, 0x18>(); // Cu, Metals Silver
      default:
                return STLFormat::colour<0x08, 0x08, 0x08>(); // Default, grey

    }
  }

  static void write_stl_vector(ostream& out, const vector3& v) {
    array<float, 3> tmp{
      static_cast<float>(v.x()), static_cast<float>(v.y()), static_cast<float>(v.z())
    };
    out.write(reinterpret_cast<const char*>(tmp.data()), tmp.size() * sizeof(float));
  }

  static void write_stl_triangle(ostream& out, const Triangle& t) {
    write_stl_vector(out, t.a);
    write_stl_vector(out, t.b);
    write_stl_vector(out, t.c);
    out.write(reinterpret_cast<const char*>(&t.col), sizeof(t.col));
  }

  static void output_stl( ostream &os, const vector<Triangle> &triangles, bool cpk_colour ) {
    static constexpr size_t StlHeaderSizeInBytes{80};
    static constexpr array<uint8_t, 10> StlGlobalColorHeader{
      'C', 'O', 'L', 'O', 'R', '=', 0xFF, 0xFF, 0xFF, 0xFF
    };
    
    // STL header
    size_t zero_bytes_in_header = StlHeaderSizeInBytes;
    if (cpk_colour) {
      os.write(reinterpret_cast<const char*>(StlGlobalColorHeader.data()), StlGlobalColorHeader.size());
      zero_bytes_in_header -= StlGlobalColorHeader.size();
    }
    for (size_t i = 0; i < zero_bytes_in_header; ++i) {
      os.put(0);
    }

    uint32_t triangle_count = triangles.size();
    os.write( reinterpret_cast<const char*>(&triangle_count), sizeof(uint32_t) );

    // go through all triangles
    for(const auto& t : triangles) {
      // write zero normal - most programs ignore the normals
      write_stl_vector(os, VZero);

      // write triangle
      write_stl_triangle(os, t);
    }

    os.flush();
  }

  bool STLFormat::WriteMolecule( OBBase* pOb, OBConversion* pConv )
  {
    OBMol* pmol = dynamic_cast< OBMol* >(pOb);
    if (pmol == nullptr) return false;

    ostream& os = *pConv->GetOutStream();

    double probe_radius = 0.;
    double scale_factor = 1.;
    uint16_t col = 0x0; // no colour colour
    bool cpk_colours = (pConv->IsOption("c") != nullptr);

    if( pConv->IsOption( "p" ) ) {
      probe_radius = atof( pConv->IsOption("p", OBConversion::OUTOPTIONS) );
      if( !isfinite(probe_radius) || probe_radius < 0. ) { probe_radius = 0.; }
    }
    if( pConv->IsOption( "s" ) ) {
      scale_factor = atof( pConv->IsOption("s", OBConversion::OUTOPTIONS) );
      if( !isfinite(scale_factor) || scale_factor < 0. ) { scale_factor = 0.; }
    }


    vector<Triangle> triangles;
    FOR_ATOMS_OF_MOL(a, *pmol) {
      const double vdwrad = scale_factor * OBElements::GetVdwRad( a->GetAtomicNum() ) + probe_radius;
      col = (cpk_colours) ? cpk_colour(  a->GetAtomicNum() ) : col;
      map_sphere( triangles, a->GetVector(), vdwrad, col );
    }

    output_stl ( os, triangles, cpk_colours );
    os.flush();
    return true;
  }

}
