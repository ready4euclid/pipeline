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
 *  @file Catalogue/Catalogue.cpp
 *
 *  @brief Methods of the class Catalogue 
 *
 *  This file contains the implementation of the methods of the class
 *  Catalogue, used to handle catalogues of astronomical sources
 *
 *  @author Federico Marulli
 *
 *  @author federico.marulli3@unbo.it
 */

#include "Catalogue.h"

using namespace cosmobl;
using namespace catalogue;


// ============================================================================

/// @cond template

template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::RandomObject>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Mock>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Halo>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Galaxy>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Cluster>);
template cosmobl::catalogue::Catalogue::Catalogue (vector<cosmobl::catalogue::Void>);

template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::RandomObject);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Mock);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Halo);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Galaxy);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Cluster);
template void cosmobl::catalogue::Catalogue::add_object (cosmobl::catalogue::Void);

template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::RandomObject>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Mock>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Halo>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Galaxy>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Cluster>);
template void cosmobl::catalogue::Catalogue::add_objects (vector<cosmobl::catalogue::Void>);

template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::RandomObject>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Mock>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Halo>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Galaxy>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Cluster>);
template void cosmobl::catalogue::Catalogue::remove_objects (vector<cosmobl::catalogue::Void>);

/// @endcond

// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const ObjType type, const vector<double> xx, const vector<double> yy, const vector<double> zz, vector<double> weight)
{
  if (weight.size() == 0)
    weight.resize(xx.size(), 1);

  for (size_t i=0; i<xx.size(); i++)
    m_sample.push_back(move(Object::Create(type, xx[i], yy[i], zz[i], weight[i])));
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const ObjType type, const vector<double> ra, const vector<double> dec, const vector<double> redshift, const Cosmology &cosm, vector<double> weight)
{
  if (weight.size() == 0)
    weight.resize(ra.size(), 1);

  for (size_t i=0; i<ra.size(); i++)
    m_sample.push_back(move(Object::Create(type, ra[i], dec[i], redshift[i], cosm, weight[i])));
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const vector<string> file, const ObjType type, const double nSub) 
{
  cout << "I'm importing the catalogue..." << endl;

  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  for (size_t dd=0; dd<file.size(); dd++) {
    string file_in = file[dd];
    cout <<"file_in = "<<file_in<<endl;
    ifstream finr (file_in.c_str()); checkIO(file_in, 1);
    
    string line; double XX, YY, ZZ;
    while (getline(finr, line)) {
      if (ran(gen)<nSub) {
	stringstream ss(line);
	ss>>XX; ss>>YY, ss>>ZZ;
        m_sample.push_back(move(Object::Create(type, XX, YY, ZZ, 1.)));
      }
    }
    
    finr.clear(); finr.close(); 
  }
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const string file_in, const Cosmology &cosm, const ObjType type, const double nSub, const double fact) 
{      
  cout << "I'm importing the catalogue: " << file_in << "..." << endl;
  
  ifstream finr(file_in.c_str());  checkIO(file_in, 1);
 
  string line;
  double RA, DEC, ZZ;

  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  while (getline(finr,line)) {
    stringstream ss(line);
    ss >> RA; ss >> DEC; ss >> ZZ;

    if (ran(gen)<nSub) 
      m_sample.push_back(move(Object::Create(type, RA, DEC, ZZ, cosm)));
  }
  
  finr.clear(); finr.close(); 
}


// ============================================================================


cosmobl::catalogue::Catalogue::Catalogue (const Catalogue catalogue, const double N_R)
{
  // ------ create the random catalogue ------ 
  
  cout << "I'm creating a random catalogue..." << endl;

  int nRandom = (int)N_R*catalogue.nObjects();

  vector<double> Lim; 
  catalogue.MinMax_var(Var::_XX_, Lim, 0);
  catalogue.MinMax_var(Var::_YY_, Lim, 0);
  catalogue.MinMax_var(Var::_ZZ_, Lim, 0);
  
  double Xmin = Lim[0]; double Xmax = Lim[1];
  double Ymin = Lim[2]; double Ymax = Lim[3];
  double Zmin = Lim[4]; double Zmax = Lim[5];

  if (Xmin>Xmax || Ymin>Ymax || Zmin>Zmax) ErrorMsg("Error in random_catalogue_box of RandomCatalogue.cpp!");

  default_random_engine gen;
  uniform_real_distribution<float> ran(0., 1.);

  double XX, YY, ZZ;

  for (int i=0; i<nRandom; i++) {
    XX = ran(gen)*(Xmax-Xmin)+Xmin;
    YY = ran(gen)*(Ymax-Ymin)+Ymin;
    ZZ = ran(gen)*(Zmax-Zmin)+Zmin;
    m_sample.push_back(move(Object::Create(_RandomObject_, XX, YY, ZZ, 1.)));
  }
}


// ============================================================================


int cosmobl::catalogue::Catalogue::Nregion () const
{ 
  vector<int> regions; 
  for (int i=0; i<nObjects(); i++) regions.push_back(m_sample[i]->region()); 
  sort(regions.begin(), regions.end());
  vector<int>::iterator it = unique(regions.begin(), regions.end());
  regions.resize(std::distance(regions.begin(), it));
  return regions.size();
}


// ============================================================================


vector<long> cosmobl::catalogue::Catalogue::get_region_list () const
{ 
  vector<long> regions; 
  for (int i=0; i<nObjects(); i++) regions.push_back(m_sample[i]->region()); 
  sort(regions.begin(), regions.end());
  vector<long>::iterator it = unique(regions.begin(), regions.end());
  regions.resize(std::distance(regions.begin(), it));
  return regions;
}


// ============================================================================


vector<double> cosmobl::catalogue::Catalogue::var (Var var_name) const
{
  vector<double> vv(m_sample.size(), 0.);
  
  switch (var_name) {

  case Var::_XX_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->xx();
    break;

  case Var::_YY_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->yy();
    break;

  case Var::_ZZ_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->zz();
    break;

  case Var::_RA_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->ra();
    break;

  case Var::_DEC_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->dec();
    break;

  case Var::_REDSHIFT_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->redshift();
    break;

  case Var::_DC_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->dc();
    break;

  case Var::_WEIGHT_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->weight();
    break;

  case Var::_MASS_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->mass();
    break;

  case Var::_MAGNITUDE_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->magnitude();
    break;

  case Var::_RICHNESS_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->richness();
    break;

  case Var::_VX_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->vx();
    break;
  
  case Var::_VY_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->vy();
    break;
  
  case Var::_VZ_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->vz();
    break;

  case Var::_REGION_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->region();
    break; 
  
  case Var::_GENERIC_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->generic();
    break;

  case Var::_RADIUS_:
    for (int i=0; i<nObjects(); i++) vv[i] = m_sample[i]->radius();
    break;

  default:
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::var of Catalogue.cpp: no such a variable in the list!");
  }
  
  return vv;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::set_var (const Var var_name, const vector<double> _var)
{
  if (m_sample.size()!=_var.size()) ErrorMsg ("Error in cosmobl::catalogue::Catalogue::set_var of Catalogue.cpp: different sizes!");
  
  switch (var_name) {

  case Var::_XX_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_xx(_var[i]);
    break;

  case Var::_YY_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_yy(_var[i]);
    break;

  case Var::_ZZ_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_zz(_var[i]);
    break;

  case Var::_RA_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_ra(_var[i]);
    break;

  case Var::_DEC_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_dec(_var[i]);
    break;

  case Var::_REDSHIFT_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_redshift(_var[i]);
    break;

  case Var::_DC_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_dc(_var[i]);
    break;

  case Var::_WEIGHT_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_weight(_var[i]);
    break;

  case Var::_MASS_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_mass(_var[i]);
    break;

  case Var::_MAGNITUDE_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_magnitude(_var[i]);
    break;

  case Var::_RICHNESS_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_richness(_var[i]);
    break;

  case Var::_VX_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_vx(_var[i]);
    break;
  
  case Var::_VY_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_vy(_var[i]);
    break;
  
  case Var::_VZ_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_vz(_var[i]);
    break;

  case Var::_REGION_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_region(_var[i]);
    break;

  case Var::_GENERIC_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_generic(_var[i]);
    break;

  case Var::_RADIUS_:
    for (int i=0; i<nObjects(); i++) m_sample[i]->set_radius(_var[i]);
    break;

  default:
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::set_var of Catalogue.cpp: no such a variable in the list!");
  }

}


// ============================================================================


void cosmobl::catalogue::Catalogue::MinMax_var (const Var var_name, vector<double> &Lim, const bool er) const
{
  if (er) Lim.erase(Lim.begin(), Lim.end());
  Lim.push_back(Min(var(var_name))); 
  Lim.push_back(Max(var(var_name)));
}


// ============================================================================


void cosmobl::catalogue::Catalogue::MinMax_var (const vector<Var> var_name, vector<vector<double> > &Lim, const bool er) const
{
  if (er) Lim.erase(Lim.begin(),Lim.end());

  for (unsigned int i=0; i<var_name.size(); i++) 
    Lim.push_back(MinMax_var(var_name[i]));
}


// ============================================================================


vector<double> cosmobl::catalogue::Catalogue::MinMax_var (const Var var_name) const
{
  vector<double> Lim; 

  Lim.push_back(Min(var(var_name))); Lim.push_back(Max(var(var_name)));

  return Lim;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::stats_var (const Var var_name, vector<double> &stats) const
{
  stats.erase(stats.begin(), stats.end());
  stats.resize(4);
  
  stats[0] = Average(var(var_name)); 
  stats[2] = Sigma(var(var_name));
  
  double first, third;
  Quartile(var(var_name), stats[1], first, third);
  
  stats[3] = third-first;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::stats_var (const vector<Var> var_name, vector<vector<double> > &stats) const
{
  stats.erase(stats.begin(),stats.end());

  for (unsigned int i=0; i<var_name.size(); i++) {

    vector<double> stats_temp;
    stats_var(var_name[i],stats_temp);

    stats.push_back(stats_temp);
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::var_distr (const Var var_name, vector<double> &_var, vector<double> &dist, const int nbin, const bool linear, const string file_out, const double Volume, const bool norm, const double V1, const double V2, const bool bin_type, const bool convolution, const double sigma) const
{ 
  distribution(_var, dist, var(var_name), var(Var::_WEIGHT_), nbin, linear, file_out, (norm) ? Volume*weightedN() : Volume, V1, V2, bin_type, convolution, sigma);
}


// ============================================================================


void cosmobl::catalogue::Catalogue::computeComovingCoordinates (const Cosmology &cosm)
{
  double red, xx, yy, zz;

  for (int i=0; i<nObjects(); i++) {

    red = redshift(i);
    m_sample[i]->set_dc(cosm.D_C(red));
    
    cartesian_coord(ra(i), dec(i), dc(i), xx, yy, zz);
    
    m_sample[i]->set_xx(xx); 
    m_sample[i]->set_yy(yy); 
    m_sample[i]->set_zz(zz);
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::computePolarCoordinates ()
{
  double ra, dec, dc;

  for (int i=0; i<nObjects(); i++) {
    polar_coord(xx(i), yy(i), zz(i), ra, dec, dc);
    m_sample[i]->set_ra(ra); 
    m_sample[i]->set_dec(dec); 
    m_sample[i]->set_dc(dc);
  }
}

// ============================================================================


void cosmobl::catalogue::Catalogue::computePolarCoordinates (const Cosmology &cosm, const double z1, const double z2)
{
  double ra, dec, dc;

  for (int i=0; i<nObjects(); i++) {
    polar_coord(xx(i), yy(i), zz(i), ra, dec, dc);
    m_sample[i]->set_ra(ra); 
    m_sample[i]->set_dec(dec); 
    m_sample[i]->set_dc(dc);
    m_sample[i]->set_redshift(cosm.Redshift(dc, z1, z2));
  }
}

// ============================================================================


void cosmobl::catalogue::Catalogue::normalizeComovingCoordinates () 
{
  for (int i=0; i<nObjects(); i++) { 
    m_sample[i]->set_xx(xx(i)/dc(i)); 
    m_sample[i]->set_yy(yy(i)/dc(i)); 
    m_sample[i]->set_zz(zz(i)/dc(i));
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::restoreComovingCoordinates ()
{
  for (int i=0; i<nObjects(); i++) {
    m_sample[i]->set_xx(xx(i)*dc(i)); 
    m_sample[i]->set_yy(yy(i)*dc(i)); 
    m_sample[i]->set_zz(zz(i)*dc(i));
  }
}


// ============================================================================


void cosmobl::catalogue::Catalogue::Order (const vector<int> vv) 
{
  int nObj = m_sample.size();

  if (int(vv.size())!=nObj) ErrorMsg("Error in cosmobl::catalogue::Catalogue::Order!");
 
  vector<shared_ptr<Object>> obj(nObj);
  
  m_index.resize(nObj);
  
  for (size_t i=0; i<vv.size(); i++) {
    m_index[i] = vv[i];
    obj[i] = m_sample[vv[i]];
  }

  m_sample = obj;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::Order () 
{ 
  int nObj = m_sample.size();
  
  vector<shared_ptr<Object>> obj(nObj);
  
  if (m_index.size() != nObj) 
    ErrorMsg("Error in Catalogue::Order() of Catalogue.cpp, order not found!");
  
  obj = m_sample;
  
  for (int i=0; i<nObj; i++) 
    m_sample[i] = obj[m_index[i]];
}


// ============================================================================


double cosmobl::catalogue::Catalogue::weightedN () const
{
  double nn = 0.;
  for (unsigned int i=0; i<m_sample.size(); i++)
    nn += m_sample[i]->weight();
  return nn;
}


// ============================================================================


void cosmobl::catalogue::Catalogue::write_comoving_coordinates (const string file_output) const
{
  ofstream fout(file_output.c_str()); checkIO(file_output, 0);
 
  for (int i=0; i<nObjects(); i++) 
    fout << xx(i) << "   " << yy(i) << "   " << zz(i) << endl;

  cout << "I wrote the file: " << file_output << endl;
  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::catalogue::Catalogue::write_obs_coordinates (const string file_output) const 
{
  ofstream fout(file_output.c_str()); checkIO(file_output, 0);

  if (!isSet(ra(0)) || !isSet(dec(0)) || !isSet(redshift(0)) || !isSet(dc(0)))
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::write_obs_coords of Catalogue.cpp! Polar coordinates are not set!");
  
  for (int i=0; i<nObjects(); i++) 
    fout << ra(i) << "   " << dec(i) << "   " << redshift(i) << "   " << dc(i) << endl;
  
  cout << "I wrote the file: " << file_output << endl;
  fout.clear(); fout.close();
}


// ============================================================================


void cosmobl::catalogue::Catalogue::write_coordinates (const string file_output) const 
{
  ofstream fout(file_output.c_str()); checkIO(file_output, 0);

  if (!isSet(ra(0)) || !isSet(dec(0)) || !isSet(redshift(0)) || !isSet(dc(0)))
    ErrorMsg("Error in cosmobl::catalogue::Catalogue::write_coordinates of Catalogue.cpp! Polar coordinates are not set!");
  
  for (int i=0; i<nObjects(); i++) 
    fout << xx(i) << "   " << yy(i) << "   " << zz(i) << "   " << ra(i) << "   " << dec(i) << "   " << redshift(i) << "   " << dc(i) << endl;
  
  cout << "I wrote the file: " << file_output << endl;
  fout.clear(); fout.close();
}


// ============================================================================


shared_ptr<Catalogue> cosmobl::catalogue::Catalogue::cut (const Var var_name, const double down, const double up, const bool excl)
{
  vector<shared_ptr<Object> > objects;
  vector<double> vvar = var(var_name);
  vector<int> w(vvar.size());

  for (size_t i=0; i<m_sample.size(); i++){
      w[i] = (excl) ? 1 : 0;
    if (vvar[i] >= down && vvar[i] < up)
      w[i] = (excl) ? 0 : 1;
  }

  for (size_t i=0; i<m_sample.size(); i++)
    if (w[i]==1)
      objects.push_back(m_sample[i]);
 
  shared_ptr<Catalogue> cat(new Catalogue{objects});
  return cat;
}


// ============================================================================


double cosmobl::catalogue::Catalogue::distance (const int i, shared_ptr<Object> obj) const
{
  return sqrt((m_sample[i]->xx()-obj->xx())*(m_sample[i]->xx()-obj->xx())+
	      (m_sample[i]->yy()-obj->yy())*(m_sample[i]->yy()-obj->yy())+
	      (m_sample[i]->zz()-obj->zz())*(m_sample[i]->zz()-obj->zz()));
}


// ============================================================================


double cosmobl::catalogue::Catalogue::angsep_xyz (const int i, shared_ptr<Object> obj) const
{ 
  return 2.*asin(0.5*sqrt((m_sample[i]->xx()-obj->xx())*(m_sample[i]->xx()-obj->xx())+
			  (m_sample[i]->yy()-obj->yy())*(m_sample[i]->yy()-obj->yy())+
			  (m_sample[i]->zz()-obj->zz())*(m_sample[i]->zz()-obj->zz())));
}
    

// ============================================================================


shared_ptr<Catalogue> cosmobl::catalogue::Catalogue::smooth (const double gridsize, const vector<Var> vars, const int SUB)
{
  shared_ptr<Catalogue> cat {new Catalogue(*this)};
  
  if (gridsize<1.e-30) return cat;
  
  double rMAX = 0.;
  
  vector<shared_ptr<Object>> sample;

  
  // ----------------------------------------------------------------------------
  // ----- subdivide the catalogue in SUB^3 sub-catalogues, to avoid ------------
  // ----- memory problems possibly caused by too large chain-mesh vectors ------
  // ----------------------------------------------------------------------------

  cout <<"Please wait, I'm subdividing the catalogue in "<<pow(SUB, 3)<<" sub-catalogues..."<<endl;
  
  vector<double> Lim;
  cat->MinMax_var(Var::_XX_, Lim, 0);
  cat->MinMax_var(Var::_YY_, Lim, 0);
  cat->MinMax_var(Var::_ZZ_, Lim, 0);

  int nx = SUB, ny = SUB, nz = SUB;
  
  double Cell_X = (Lim[1]-Lim[0])/nx;
  double Cell_Y = (Lim[3]-Lim[2])/ny;
  double Cell_Z = (Lim[5]-Lim[4])/nz;

  for (int i=0; i<cat->nObjects(); i++) {
    int i1 = min(int((cat->xx(i)-Lim[0])/Cell_X), nx-1);
    int j1 = min(int((cat->yy(i)-Lim[2])/Cell_Y), ny-1);
    int z1 = min(int((cat->zz(i)-Lim[4])/Cell_Z), nz-1);
    int index = z1+nz*(j1+ny*i1);
    cat->object(i)->set_region(index);
  }
  
  vector<long> region_list = cat->get_region_list();
  int nRegions = region_list.size();
  
  vector<shared_ptr<Catalogue>> subSamples(nRegions);
  
  for (int i=0; i<nRegions; i++) {
    double start = region_list[i];
    double stop = start+1;
    subSamples[i] = cut(Var::_REGION_, start, stop);
  }

  
  // -----------------------------------------
  // ----- smooth all the sub-catalogues -----
  // -----------------------------------------

  for (int rr=0; rr<nRegions; rr++) {
    
    vector<double> _xx = subSamples[rr]->var(Var::_XX_), _yy = subSamples[rr]->var(Var::_YY_), _zz = subSamples[rr]->var(Var::_ZZ_);
    
    ChainMesh3D ll(gridsize, _xx, _yy, _zz, rMAX, (long)-1.e5, (long)1.e5);
   
    
    for (int i=0; i<ll.nCell(); i++) {
      vector<long> list = ll.get_list(i);
     
      int nObj = list.size();
      if (nObj>0) {
	
	double XX = 0., YY = 0., ZZ = 0., RA = 0., DEC = 0., REDSHIFT = 0., WEIGHT = 0.;

	for (size_t j=0; j<list.size(); j++) {
	  XX += subSamples[rr]->xx(list[j]);
	  YY += subSamples[rr]->yy(list[j]);
	  ZZ += subSamples[rr]->zz(list[j]);
	  RA += subSamples[rr]->ra(list[j]);
	  DEC += subSamples[rr]->dec(list[j]);
	  REDSHIFT += subSamples[rr]->redshift(list[j]);
	  WEIGHT += subSamples[rr]->weight(list[j]);
	}
	
	shared_ptr<Object> obj{new Object()};
	XX /= nObj; obj->set_xx(XX);
	YY /= nObj; obj->set_yy(YY);
	ZZ /= nObj; obj->set_zz(ZZ);
	RA /= nObj; obj->set_ra(RA);
	DEC /= nObj; obj->set_dec(DEC);
	REDSHIFT /= nObj; obj->set_redshift(REDSHIFT);
	obj->set_weight(WEIGHT);
	sample.push_back(obj);

      }
    }
    
  }
  
  shared_ptr<Catalogue> cat_new(new Catalogue{sample});
  return cat_new;
}


// ============================================================================


int cosmobl::catalogue::Catalogue::nObjects_condition (const Var var_name, const double down, const double up, const bool excl)
{
  int nObjw = 0;
  vector<double> vvar = var(var_name);

  for (size_t i=0; i<m_sample.size(); i++)

    if (vvar[i] >= down && vvar[i] < up) 
      nObjw ++;
    
  nObjw = (excl) ? weightedN()-nObjw : nObjw;

  return nObjw;
}


// ============================================================================


double cosmobl::catalogue::Catalogue::weightedN_condition (const Var var_name, const double down, const double up, const bool excl)
{
  double nObjw = 0;
  vector<double> vvar = var(var_name);

  for (size_t i=0; i<m_sample.size(); i++)
    if (vvar[i] >= down && vvar[i] < up)
      nObjw += weight(i);
    
  nObjw = (excl) ? weightedN()-nObjw : nObjw;

  return nObjw;
}
