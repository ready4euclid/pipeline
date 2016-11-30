/********************************************************************
 *  Copyright (C) 2015 by Federico Marulli and Alfonso Veropalumbo  *
 *  federico.marulli3@unibo.it                                      *
 *                                                                  *
 *  This program is free software; you can redistribute it and/or   *
 *  modify it under the terms of the GNU General Public License as  *
 *  published by the Free Software Foundation; either version 2 of  *
 *  the License, or (at your option) any later version.             *
 *                                                                  *
 *  This program is distributed in the hope that it will be useful, *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of  *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the   *
 *  GNU General Public License for more details.                    *
 *                                                                  *
 *  You should have received a copy of the GNU General Public       *
 *  License along with this program; if not, write to the Free      *
 *  Software Foundation, Inc.,                                      *
 *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.       *
 ********************************************************************/

/**
 *  @file Headers/Objects/Mock.h
 *
 *  @brief The class Mock 
 *
 *  This file defines the interface of the class Mock, used to handle
 *  objects of type <EM> mock object </EM>
 *
 *  @authors Federico Marulli, Alfonso Veropalumbo
 *
 *  @authors federico.marulli3@unbo.it, alfonso.veropalumbo@unibo.it
 */

#ifndef __MOCK__
#define __MOCK__ 


// ===================================================================================================


namespace cosmobl {

  namespace catalogue {
    
    /**
     *  @class Mock Mock.h "Headers/Lib/Mock.h"
     *
     *  @brief The class Mock
     *
     *  This class is used to handle objects of type <EM> mock object
     *  </EM>
     */
    class Mock : public Halo {

    private :

      /// generic variable of the mock object
      double m_generic; 

    
    public:

      /**
       *  @name Constructors/destructors
       */
      ///@{

      /**
       *  @brief default constructor
       *  @return object of class Mock
       */
      Mock ()
	: Halo(), m_generic(par::defaultDouble) {}

      /**
       *  @brief constructor that uses comoving coordinates
       *
       *  @param coord structure containing the comoving coordinates
       *  {x, y, z}
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param generic the mock generic variable
       *
       *  @return object of class Mock
       */
      Mock (const comovingCoordinates coord, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double generic=par::defaultDouble) 
	: Halo(coord, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_generic(generic) {}

      /**
       *  @brief constructor that uses comoving coordinates and a
       *  cosmological model to estimate the redshift
       *
       *  @param coord structure containing the comoving coordinates
       *  {x, y, z}
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param z1_guess minimum prior on the redshift
       *
       *  @param z2_guess maximum prior on the redshift 
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *   
       *  @param generic the mock generic variable
       *
       *  @return object of class Mock
       */
      Mock (const comovingCoordinates coord, const cosmology::Cosmology &cosm, const double z1_guess=0., const double z2_guess=10., const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double generic=par::defaultDouble) 
	: Halo(coord, cosm, z1_guess, z2_guess, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_generic(generic) {}

      /**
       *  @brief constructor that uses observed coordinates in radians
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param generic the mock generic variable
       *
       *  @return object of class Mock
       */
      Mock (const observedCoordinates coord, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double generic=par::defaultDouble) 
	: Halo(coord, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_generic(generic) {}
      
      /**
       *  @brief constructor that uses observed coordinates in any
       *  angular units
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param generic the mock generic variable
       *
       *  @return object of class Mock
       */
      Mock (const observedCoordinates coord, const CoordUnits inputUnits, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double generic=par::defaultDouble) 
	: Halo(coord, inputUnits, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_generic(generic) {}
      
      /**
       *  @brief constructor that uses observed coordinates in radians
       *  and a cosmological model to estimate the comoving
       *  coordinates
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshitf}
       *
       *  @param cosm object of class Cosmology, used to estimate
       *  comoving distances
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param generic the mock generic variable
       *
       *  @return object of class Mock
       */
      Mock (const observedCoordinates coord, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double generic=par::defaultDouble) 
	: Halo(coord, cosm, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_generic(generic) {}

      /**
       *  @brief constructor that uses observed coordinates and a
       *  cosmological model to estimate the comoving coordinates
       *
       *  @param coord structure containing the observed coordinates
       *  {R.A., Dec, redshift}
       *
       *  @param inputUnits the units of the input coordinates
       *
       *  @param cosm object of class Cosmology, used to estimate comoving distances
       *
       *  @param weight weight
       *
       *  @param region region, used e.g. for jackknife and bootstrap
       *
       *  @param field the field where the object has been observed
       *
       *  @param x_displacement the displacement along the x-axis
       *
       *  @param y_displacement the displacement along the y-axis
       *
       *  @param z_displacement the displacement along the z-axis
       *
       *  @param vx halo peculiar velocity along the x direction
       *
       *  @param vy halo peculiar velocity along the y direction
       *
       *  @param vz halo peculiar velocity along the z direction
       *
       *  @param mass the mock mass
       *
       *  @param generic the mock generic variable
       *
       *  @return object of class Mock
       */
      Mock (const observedCoordinates coord, const CoordUnits inputUnits, const cosmology::Cosmology &cosm, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double generic=par::defaultDouble) 
	: Halo(coord, inputUnits, cosm, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_generic(generic) {}

      /**
       *  @brief constructor that uses both comoving and observed coordinates
       *  @param xx comoving coordinate
       *  @param yy comoving coordinate
       *  @param zz comoving coordinate 
       *  @param ra Right Ascension
       *  @param dec Declination
       *  @param redshift redshift
       *  @param weight weight   
       *  @param region region, used e.g. for jackknife and bootstrap
       *  @param field the field where the object has been observed
       *  @param x_displacement the displacement along the x-axis
       *  @param y_displacement the displacement along the y-axis
       *  @param z_displacement the displacement along the z-axis
       *  @param vx halo peculiar velocity along the x direction
       *  @param vy halo peculiar velocity along the y direction
       *  @param vz halo peculiar velocity along the z direction
       *  @param mass the mock mass
       *  @param generic the mock generic variable
       *
       *  @return object of class Mock
       */
      Mock (const double xx, const double yy, const double zz, const double ra, const double dec, const double redshift, const double weight=1., const long region=par::defaultLong, const string field=par::defaultString, const double x_displacement=par::defaultDouble, const double y_displacement=par::defaultDouble, const double z_displacement=par::defaultDouble, const double vx=par::defaultDouble, const double vy=par::defaultDouble, const double vz=par::defaultDouble, const double mass=par::defaultDouble, const double generic=par::defaultDouble) 
	: Halo(xx, yy, zz, ra, dec, redshift, weight, region, field, x_displacement, y_displacement, z_displacement, vx, vy, vz, mass), m_generic(generic) {}
      
      /**
       *  @brief default destructor
       *  @return none
       */
      ~Mock () = default;

      ///@}
  
  
      /**
       *  @name Member functions used to get the private members 
       */
      ///@{
    
      /**
       *  @brief get the private member Mock::m_generic
       *  @return the generic variable of the mock object
       */
      double generic () const override { return m_generic; }

      ///@}
  
  
      /**
       *  @name Member functions used to set the private members 
       */
      ///@{
    
      /**
       *  @brief set the private member Mock::m_generic
       *  @param generic the generic variable of the mock object
       *  @return none
       */
      void set_generic (const double generic=par::defaultDouble) override { m_generic = generic; }
    
      ///@}

    };
  }
}

#endif
