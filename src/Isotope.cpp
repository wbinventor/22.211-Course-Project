/*
 * Isotope.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Isotope.h"


/**
 * Isotope constructor sets default values for some isotope properties
 */
Isotope::Isotope() {

	/* Default atomic number and number densities and temperature */
	_isotope_name = (char*)"";
	_A = 1;
	_N = 1;
	_T = 300;
	_mu_avg = 2.0 / (3.0 * _A);
	_kB = 8.617332E-5;             /* boltzmann's constant (ev / K) */
	_alpha = float(_A-1)/float(_A+1) * float(_A-1)/float(_A+1);
	_eta = (float(_A)+1.0) / (2.0 * sqrt(float(_A)));
	_rho = (float(_A)-1.0) / (2.0 * sqrt(float(_A)));

	/* By default this isotope has no cross-sections */
	_num_capture_xs = 0;
	_num_scatter_xs = 0;
	_num_elastic_xs = 0;
	_num_inelastic_xs = 0;
	_num_fission_xs = 0;
	_num_total_xs = 0;

	/* By default the thermal scattering cdfs have not been initialized */
	_num_thermal_cdfs = 0;
	_num_thermal_cdf_bins = 0;
}


/**
 * Isotope destructor deletes array of cross section values that
 * have been assigned to this isotope
 */
Isotope::~Isotope() {
	if (_num_capture_xs != 0) {
		delete [] _capture_xs;
		delete [] _capture_xs_energies;
	}
	if (_num_scatter_xs != 0) {
		delete [] _scatter_xs;
		delete [] _scatter_xs_energies;
	}
	if (_num_elastic_xs != 0) {
		delete [] _elastic_xs;
		delete [] _elastic_xs_energies;
	}
	if (_num_inelastic_xs != 0) {
		delete [] _inelastic_xs;
		delete [] _inelastic_xs_energies;
	}
	if (_num_fission_xs != 0) {
		delete [] _fission_xs;
		delete [] _fission_xs_energies;
	}
	if (_num_total_xs != 0) {
		delete [] _total_xs;
		delete [] _total_xs_energies;
	}
	if (_num_thermal_cdfs != 0) {
		delete [] _thermal_dist;
		for (int i=0; i < _num_thermal_cdfs; i++)
			delete [] _thermal_cdfs[i];
		delete [] _thermal_cdfs;
		delete [] _E_to_kT;
		delete [] _Eprime_to_E;
	}
}


/**
 * Returns the name of the of isotope
 * @return character array with name of isotope
 */
char* Isotope::getIsotopeType() const {
	return _isotope_name;
}


/**
 * Returns the atomic number of this isotope
 * @return the atomic number
 */
int Isotope::getA() const {
    return _A;
}


/**
 * Returns the alpha ((A-1)/(A+1))^2 values for this isotope
 * @return alpha
 */
float Isotope::getAlpha() const {
    return _alpha;
}


/**
 * Returns the number density for this isotope
 * @return the number density
 */
float Isotope::getN() const {
    return _N;
}


/**
 * Return the temperature (Kelvin) for this isotope
 * @return the temperature of this isotope
 */
float Isotope::getTemperature() const {
	return _T;
}


/**
 * Return the average value of the cosine of theta for this isotope
 * in a scattering collision
 * @return the average for mu
 */
float Isotope::getMuAverage() const {
	return _mu_avg;
}


/**
 * Returns a capture cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the capture cross-section (barns)
 */
float Isotope::getCaptureXS(float energy) const{

	if (_num_capture_xs == 0)
		return 0.0;

	/* Use linear interpolation to find the capture cross-section */
	return linearInterp(_capture_xs_energies, _capture_xs,
									_num_capture_xs, energy);
}


/**
 * Returns a capture cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the capture cross-section (barns)
 */
float Isotope::getCaptureXS(int energy_index) const {

	if (_num_capture_xs == 0)
		return 0.0;

	else if (energy_index > _num_capture_xs)
		log_printf(ERROR, "Unable to retrieve capture xs for"
				" isotope %s since the energy index %d is out of"
				" bounds", _isotope_name, energy_index);

	return _capture_xs[energy_index];
}


/**
 * Returns a one group scattering cross-section value
 * @return the scattering cross-section (barns)
 */
float Isotope::getOneGroupCaptureXS() const {
	if (_num_capture_xs == 0)
		return 0.0;

	return _capture_xs[0];
}


/**
 * Returns an scattering cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the scattering cross-section (barns)
 */
float Isotope::getScatterXS(float energy) const {
	if (_num_scatter_xs == 0 && _num_inelastic_xs == 0 && _num_elastic_xs == 0)
		return 0.0;

	/* If this isotope has both inelastic and elastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs != 0 &&
												_num_elastic_xs != 0)
		return (getInelasticXS(energy) + getElasticXS(energy));

	/* If this isotope only has elastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs == 0 &&
												_num_elastic_xs != 0)
		return getElasticXS(energy);

	/* If this isotope only has inelastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs != 0 &&
												_num_elastic_xs == 0)
		return getInelasticXS(energy);

	/* Use linear interpolation to find the scattering cross-section */
	return linearInterp(_scatter_xs_energies, _scatter_xs,
									_num_scatter_xs, energy);
}


/**
 * Returns a scatter cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the scatter cross-section (barns)
 */
float Isotope::getScatterXS(int energy_index) const {

	if (_num_scatter_xs == 0 && _num_inelastic_xs == 0 && _num_elastic_xs == 0)
		return 0.0;

	/* If this isotope has both inelastic and elastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs != 0 &&
												_num_elastic_xs != 0)
		return (getInelasticXS(energy_index) + getElasticXS(energy_index));

	/* If this isotope only has elastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs == 0 &&
												_num_elastic_xs != 0)
		return getElasticXS(energy_index);

	/* If this isotope only has inelastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs != 0 &&
												_num_elastic_xs == 0)
		return getInelasticXS(energy_index);

	else if (energy_index > _num_scatter_xs)
		log_printf(ERROR, "Unable to retrieve scatter xs for"
				" isotope %s since the energy index %d is out of"
				" bounds", _isotope_name, energy_index);

	return _scatter_xs[energy_index];
}


/**
 * Returns a one group scattering cross-section value
 * @return the scattering cross-section (barns)
 */
float Isotope::getOneGroupScatterXS() const {

	if (_num_scatter_xs == 0 && _num_inelastic_xs == 0 && _num_elastic_xs == 0)
		return 0.0;

	/* If this isotope has both inelastic and elastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs != 0
											&& _num_elastic_xs != 0)
		return (getOneGroupInelasticXS() + getOneGroupElasticXS());

	/* If this isotope only has elastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs == 0
											&& _num_elastic_xs != 0)
		return getOneGroupElasticXS();

	/* If this isotope only has inelastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs != 0
											&& _num_elastic_xs == 0)
		return getOneGroupInelasticXS();

	/* Otherwise, return the total scatter xs */
	return _scatter_xs[0];
}

/**
 * Returns the type of angular scattering distribution for this isotope
 * (ISOTROPIC_CM or ISOTROPIC_LAB)
 * @return the type of angular scattering distribution
 */
scatterAngleType Isotope::getScatterAngleType() const {

	if (_num_scatter_xs == 0 && _num_inelastic_xs == 0 && _num_elastic_xs == 0)
		log_printf(ERROR, "Cannot return a scatter angle type"
				"for isotope %s since it has not been set", _isotope_name);

	/* If this isotope has both inelastic and elastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs != 0
											&& _num_elastic_xs != 0)
		return getInelasticAngleType();

	/* If this isotope only has elastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs == 0
											&& _num_elastic_xs != 0)
		return getElasticAngleType();

	/* If this isotope only has inelastic scattering */
	else if (_num_scatter_xs == 0 && _num_inelastic_xs != 0
											&& _num_elastic_xs == 0)
		return getInelasticAngleType();;

	/* Otherwise, return the total scatter angle type */
	return _scatter_angle;
}


/**
 * Returns an inelastic scattering cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the inelastic scattering cross-section (barns)
 */
float Isotope::getInelasticXS(float energy) const {

	if (_num_inelastic_xs == 0)
		return 0.0;

	/* Use linear interpolation to find the inelastic scatter cross-section */
	return linearInterp(_inelastic_xs_energies, _inelastic_xs,
									_num_inelastic_xs, energy);
}


/**
 * Returns an inelastic cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the inelastic cross-section (barns)
 */
float Isotope::getInelasticXS(int energy_index) const {

	if (_num_inelastic_xs == 0)
		return 0.0;

	else if (energy_index > _num_inelastic_xs)
		log_printf(ERROR, "Unable to retrieve inelastic xs for"
				" isotope %s since the energy index %d is out of"
				" bounds", _isotope_name, energy_index);

	return _inelastic_xs[energy_index];
}


/**
 * Returns a one group inelastic scattering cross-section value
 * @return the inelastic cross-section (barns)
 */
float Isotope::getOneGroupInelasticXS() const {

	if (_num_inelastic_xs == 0)
		return 0.0;

	return _inelastic_xs[0];

}


/**
 * Returns the type of angular inelastic scattering distribution for
 * this isotope (ISOTROPIC_CM or ISOTROPIC_LAB)
 * @return the type of angular scattering distribution
 */
scatterAngleType Isotope::getInelasticAngleType() const {
	if (_num_inelastic_xs == 0)
		log_printf(ERROR, "Cannot return an inelastic angle type"
				"for isotope %s since it has not been set", _isotope_name);

	return _inelastic_angle;
}


/**
 * Returns an elastic scattering cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the elastic scattering cross-section (barns)
 */
float Isotope::getElasticXS(float energy) const {

	if (_num_elastic_xs == 0)
		return 0.0;

	/* Use linear interpolation to find the elastic scatter cross-section */
	return linearInterp(_elastic_xs_energies, _elastic_xs,
									_num_elastic_xs, energy);
}


/**
 * Returns an elastic cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the elastic cross-section (barns)
 */
float Isotope::getElasticXS(int energy_index) const {

	if (_num_elastic_xs == 0)
		return 0.0;

	else if (energy_index > _num_elastic_xs)
		log_printf(ERROR, "Unable to retrieve elastic xs for"
				" isotope %s since the energy index %d is out of"
				" bounds", _isotope_name, energy_index);

	return _elastic_xs[energy_index];
}


/**
 * Returns a one group elastic scattering cross-section value
 * @return the elastic cross-section (barns)
 */
float Isotope::getOneGroupElasticXS() const {

	if (_num_elastic_xs == 0)
		return 0.0;

	return _elastic_xs[0];
}


/**
 * Returns the type of angular elastic scattering distribution for
 * this isotope (ISOTROPIC_CM or ISOTROPIC_LAB)
 * @return the type of angular scattering distribution
 */
scatterAngleType Isotope::getElasticAngleType() const {

	if (_num_elastic_xs == 0)
		log_printf(ERROR, "Cannot return an elastic angle type"
				"for isotope %s since it has not been set", _isotope_name);

	return _elastic_angle;
}


/**
 * Returns a fission cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the fission cross-section (barns)
 */
float Isotope::getFissionXS(float energy) const{

	if (_num_fission_xs == 0)
		return 0.0;

	/* Use linear interpolation to find the fission cross-section */
	return linearInterp(_fission_xs_energies, _fission_xs,
									_num_fission_xs, energy);
}


/**
 * Returns a fission cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the fission cross-section (barns)
 */
float Isotope::getFissionXS(int energy_index) const {

	if (_num_fission_xs == 0)
		return 0.0;

	else if (energy_index > _num_fission_xs)
		log_printf(ERROR, "Unable to retrieve fission xs for"
				" isotope %s since the energy index %d is out of"
				" bounds", _isotope_name, energy_index);

	return _fission_xs[energy_index];
}


/**
 * Returns a one group fission cross-section value
 * @return the fission cross-section (barns)
 */
float Isotope::getOneGroupFissionXS() const {
	if (_num_fission_xs == 0)
		return 0.0;

	return _fission_xs[0];
}


/**
 * Returns an absorption (capture plus fission) cross-section
 * value for a certain energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the absorption cross-section (barns)
 */
float Isotope::getAbsorbXS(float energy) const {
	return (getFissionXS(energy) + getCaptureXS(energy));
}


/**
 * Returns an absorption cross-section value for a certain index into
 * the energy. This is meant to be used by a Material with rescaled
 * cross-sections
 * @param energy_index the index into the energy array
 * @return the absorption cross-section (barns)
 */
float Isotope::getAbsorbXS(int energy_index) const {
	return (getCaptureXS(energy_index) + getFissionXS(energy_index));
}


/**
 * Returns a one group absorption (capture plus fission) cross-section value
 * @return the absorption cross-section (barns)
 */
float Isotope::getOneGroupAbsorbXS() const {
	return (getOneGroupFissionXS() + getOneGroupCaptureXS());
}


/**
 * Returns a total scattering cross-section value for a certain
 * energy based on linear interpolation
 * @param energy the energy (eV) of interest
 * @return the total cross-section (barns)
 */
float Isotope::getTotalXS(float energy) const {

	/* If the total xs has been defined explicitly, use it to
	 * linearly interpolate to find the total cross-section */
	if (_num_total_xs != 0)
		return linearInterp(_total_xs_energies, _total_xs,
									_num_total_xs, energy);

	/* Otherwise loop over all xs which have been defined and
	 * add them to a total xs */
	else {

		float total_xs = 0;

		std::map<collisionType, float(Isotope::*)(float)
										const>::const_iterator iter;
		for (iter = _xs_handles.begin(); iter!= _xs_handles.end(); ++iter)
			total_xs += (this->*iter->second)(energy);

		return total_xs;
	}
}


/**
 * Returns a total scattering cross-section value for a certain index
 * into the energy array
 * @param energy the energy (eV) of interest
 * @return the total cross-section (barns)
 */
float Isotope::getTotalXS(int energy_index) const {

	/* If the total xs has been defined explicitly, use it */
	if (_num_total_xs != 0) {

		if (energy_index > _num_total_xs)
			log_printf(ERROR, "Unable to retrieve total xs for"
					" isotope %s since the energy index %d is out of"
					" bounds", _isotope_name, energy_index);

		return _total_xs[energy_index];
	}

	/* Otherwise loop over all xs which have been defined and
	 * add them to a total xs */
	else {
		float total_xs = 0;
		total_xs += getCaptureXS(energy_index);
		total_xs += getScatterXS(energy_index);
		total_xs += getFissionXS(energy_index);
		return total_xs;
	}
}


/**
 * Returns a one group total scattering cross-section value
 * @return the elastic cross-section (barns)
 */
float Isotope::getOneGroupTotalXS() const {

	/* If the total xs has been defined explicitly, use it */
	if (_num_total_xs != 0)
		return _total_xs[0];

	/* Otherwise loop over all xs which have been defined and
	 * add them to a total xs */
	else {

		float total_xs = 0;

		std::map<collisionType, float(Isotope::*)(float)
										const>::const_iterator iter;
		for (iter = _xs_handles.begin(); iter != _xs_handles.end(); ++iter)
			total_xs += (this->*iter->second)(1.0);

		return total_xs;
	}
}


/**
 * Returns the transport cross-section for this isotope at a particular energy
 * @param energy the energy (eV) of interest
 * @return the transport cross-section (barns)
 */
float Isotope::getTransportXS(float energy) const {
	return (getTotalXS(energy) - _mu_avg * getScatterXS(energy));
}


/**
 * Returns the transport cross-section for this isotope at a particular index
 * into its energy array
 * @param energy_index the index into the energy array
 * @return the transport cross-section (barns)
 */
float Isotope::getTransportXS(int energy_index) const {
	return (getTotalXS(energy_index) - _mu_avg * getScatterXS(energy_index));
}


/**
 * Computes the scattering energy from an ineleastic scattering collision
 * at a certain energy
 * @param energy of interest (eV)
 * @return the sampled scattered eprime for that energy
 */
float Isotope::getInelasticScatterEnergy(float energy) {

	if (_num_inelastic_xs == 0)
		log_printf(ERROR, "Cannot return an inelastic cross-section"
				"for isotope %s since it has not been set", _isotope_name);

	std::map<std::pair<float, float>, std::pair<std::pair<float*, float*>,
													int> >::iterator iter;
	std::pair<std::pair<float*, float*>, int> cdf;
	std::pair<float, float> cdf_bounds;
	bool found_cdf = false;

	/* Find the appropriate cdf for this energy */
	for (iter = _inelastic_cdfs.begin(); iter !=_inelastic_cdfs.end(); ++iter){

		cdf_bounds = iter->first;

		if (energy >= cdf_bounds.first && energy <= cdf_bounds.second) {
			cdf = iter->second;
			found_cdf = true;
			break;
		}
	}

	/* If no appropriate cdf was found, return the incoming energy */
	if (!found_cdf)
		return energy;

	/* If cdf was found, interpolate to find the outgoing energy */
	/* First find the cdf prob value for the incoming energy and normalize our
	 * random number test value by it */
	float prob = linearInterp(cdf.first.first, cdf.first.second, cdf.second,
																	energy);
	float test = (float(rand()) / RAND_MAX) * prob;

	/* Interpolate to find eprime based on our normalized random number */
	float energy_prime = linearInterp(cdf.first.second, cdf.first.first,
															cdf.second, test);

	return energy_prime;
}


/**
 * This method returns true if the thermal scattering distributions
 * for this isotope have been initialized, and false otherwise
 * @return boolean if the thermal scattering distributions exist
 */
bool Isotope::usesThermalScattering() {

	if (_num_thermal_cdfs == 0)
		return false;
	else
		return true;
}


/**
 * Set the isotope name
 * @param istope a character array of the isotopes name
 */
void Isotope::setIsotopeType(char* isotope) {
	_isotope_name = isotope;
}


/**
 * Set the atomic number and update alpha, eta and rho
 * @param A atomic number
 */
void Isotope::setA(int A) {
    _A = A;
	_alpha = float(_A-1)/float(_A+1) * float(_A-1)/float(_A+1);
	_eta = (float(_A)+1.0) / (2.0 * sqrt(float(_A)));
	_rho = (float(_A)-1.0) / (2.0 * sqrt(float(_A)));
	_mu_avg = 2.0 / (3.0 * _A);
}


/**
 * Set the number density (at/cm^3)
 * @param N number density (at/cm^3)
 */
void Isotope::setN(float N) {
    _N = N;
}


/* Set the temperature (Kelvin)
 * @param T the temperature (Kelvin)
 */
void Isotope::setTemperature(float T) {
	_T = T;
}


/**
 * Load a cross-section from an ASCII file into this isotope
 * @param filename the file with the cross-section values
 * @param type the type of cross-section
 * @param angle_type the type of angle (only used for scattering)
 * @param delimiter the character between data values in file
 */
void Isotope::loadXS(char* filename, collisionType type, char* delimiter) {

	/* Find the number of cross-section values in the file */
	int num_xs_values = getNumCrossSectionDataPoints(filename);

	/* Initialize data structures to store cross-section values */
	float* energies = new float[num_xs_values];
	float* xs_values = new float[num_xs_values];

	/* Parse the file into the data structures */
	parseCrossSections(filename, energies, xs_values, num_xs_values,
														delimiter);

	/* Set this isotope's appropriate cross-section using the data
	 * structures */
	if (type == CAPTURE)
		setCaptureXS(xs_values, energies, num_xs_values);
	else if (type == SCATTER)
		setScatterXS(xs_values, energies, num_xs_values, ISOTROPIC_LAB);
	else if (type == INELASTIC)
		setInelasticXS(xs_values, energies, num_xs_values, ISOTROPIC_LAB);
	else if (type == ELASTIC)
		setElasticXS(xs_values, energies, num_xs_values, ISOTROPIC_LAB);
	else if (type == FISSION)
		setFissionXS(xs_values, energies, num_xs_values);
	else
		setTotalXS(xs_values, energies, num_xs_values);

	return;
}



/**
 * Set the capture cross-section for this isotope
 * @param capture_xs a float array of microscopic capture xs (barns)
 * @param capture_xs_energies a float array of energies (eV)
 * @param num_capture_xs the number of capture xs values
 */
void Isotope::setCaptureXS(float* capture_xs, float* capture_xs_energies,
													int num_capture_xs) {
    _capture_xs = capture_xs;
    _capture_xs_energies = capture_xs_energies;
    _num_capture_xs = num_capture_xs;
    float (Isotope::*func)(float) const;
    func = &Isotope::getCaptureXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    											const>(CAPTURE, func));
}


/**
 * Set the one group capture cross-section for this isotope
 * @param capture_xs the one group capture cross-section (barns)
 */
void Isotope::setOneGroupCaptureXS(float capture_xs) {
	_capture_xs = new float[1];
	_capture_xs_energies = new float[1];
	_capture_xs[0] = capture_xs;
	_capture_xs_energies[0] = 1.0;
	_num_capture_xs = 1;
    float (Isotope::*func)(float) const;
    func = &Isotope::getCaptureXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    											const>(CAPTURE, func));
}


/**
 * Set the scattering cross-section for this isotope
 * @param scatter_xs a float array of microscopic scattering xs (barns)
 * @param scatter_xs_energies a float array of energies (eV)
 * @param num_scatter_xs the number of scattering xs values
 * @param type the type of angular scattering distribution
 */
void Isotope::setScatterXS(float* scatter_xs, float* scatter_xs_energies,
							int num_scatter_xs, scatterAngleType type) {
    _scatter_xs = scatter_xs;
    _scatter_xs_energies = scatter_xs_energies;
    _num_scatter_xs = num_scatter_xs;
    _scatter_angle = type;
    float (Isotope::*func)(float) const;
    func = &Isotope::getScatterXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    												const>(SCATTER, func));
}


/**
 * Set the one group scattering cross-section for this isotope
 * @param scatter_xs the one group scattering cross-section (barns)
 */
void Isotope::setOneGroupScatterXS(float scatter_xs, scatterAngleType type) {
	_scatter_xs = new float[1];
	_scatter_xs_energies = new float[1];
	_scatter_xs[0] = scatter_xs;
	_scatter_xs_energies[0] = 1.0;
	_num_scatter_xs = 1;
	_scatter_angle = type;
    float (Isotope::*func)(float) const;
    func = &Isotope::getScatterXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    												const>(SCATTER, func));
}


/**
 * Set the type of angular scattering distribution for this isotope
 * @param type the angular scattering distribution type
 */
void Isotope::setScatterAngleType(scatterAngleType type) {
	_scatter_angle = type;
}


/**
 * Set the inelastic cross-section for this isotope
 * @param inelastic_xs a float array of microscopic inelastic xs (barns)
 * @param inelastic_xs_energies a float array of energies (eV)
 * @param num_inelastic_xs the number of inelastic xs values
 * @param type the type of angular scattering distribution
 */
void Isotope::setInelasticXS(float* inelastic_xs,
					float* inelastic_xs_energies, int num_inelastic_xs,
												scatterAngleType type) {

	_inelastic_xs = inelastic_xs;
    _inelastic_xs_energies = inelastic_xs_energies;
    _num_inelastic_xs = num_inelastic_xs;
    _inelastic_angle = type;
    float (Isotope::*func)(float) const;
    func = &Isotope::getInelasticXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    												const>(INELASTIC, func));
}


/**
 * Set the one group inelastic scattering cross-section for this isotope
 * @param inelastic_xs the one group inelastic scattering cross-section (barns)
 */
void Isotope::setOneGroupInelasticXS(float inelastic_xs,
													scatterAngleType type) {
	_inelastic_xs = new float[1];
	_inelastic_xs_energies = new float[1];
	_inelastic_xs[0] = inelastic_xs;
	_inelastic_xs_energies[0] = 1.0;
	_num_inelastic_xs = 1;
	_inelastic_angle = type;
    float (Isotope::*func)(float) const;
    func = &Isotope::getInelasticXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    											const>(INELASTIC, func));
}


/**
 * Set the type of angular inelastic scattering distribution for this isotope
 * @param type the angular inelastic scattering distribution type
 */
void Isotope::setInelasticAngleType(scatterAngleType type) {
	_inelastic_angle = type;
}


/**
 * Set the elastic cross-section for this isotope
 * @param elastic_xs a float array of microscopic elastic xs (barns)
 * @param elastic_xs_energies a float array of energies (eV)
 * @param num_elastic_xs the number of elastic xs values
 * @param type the type of angular scattering distribution
 */
void Isotope::setElasticXS(float* elastic_xs, float* elastic_xs_energies,
								int num_elastic_xs, scatterAngleType type) {
    _elastic_xs = elastic_xs;
    _elastic_xs_energies = elastic_xs_energies;
    _num_elastic_xs = num_elastic_xs;
    _elastic_angle = type;
    float (Isotope::*func)(float) const;
    func = &Isotope::getElasticXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    												const>(ELASTIC, func));
}


/**
 * Set the one group elastic scattering cross-section for this isotope
 * @param elastic_xs the one group elastic scattering cross-section (barns)
 */
void Isotope::setOneGroupElasticXS(float elastic_xs, scatterAngleType type) {
	_elastic_xs = new float[1];
	_elastic_xs_energies = new float[1];
	_elastic_xs[0] = elastic_xs;
	_elastic_xs_energies[0] = 1.0;
	_num_elastic_xs = 1;
	_elastic_angle = type;
    float (Isotope::*func)(float) const;
    func = &Isotope::getElasticXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    												const>(ELASTIC, func));
}


/**
 * Set the type of angular elastic scattering distribution for this isotope
 * @param type the angular elastic scattering distribution type
 */
void Isotope::setElasticAngleType(scatterAngleType type) {
	_elastic_angle = type;
}


/**
 * Set the fission cross-section for this isotope
 * @param fission_xs a float array of microscopic fission xs (barns)
 * @param fission_xs_energies a float array of energies (eV)
 * @param num_fission_xs the number of fission xs values
 */
void Isotope::setFissionXS(float* fission_xs, float* fission_xs_energies,
													int num_fission_xs) {
    _fission_xs = fission_xs;
    _fission_xs_energies = fission_xs_energies;
    _num_fission_xs = num_fission_xs;
    float (Isotope::*func)(float) const;
    func = &Isotope::getFissionXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    											const>(FISSION, func));
}


/**
 * Set the one group fission cross-section for this isotope
 * @param fission_xs the one group fission cross-section (barns)
 */
void Isotope::setOneGroupFissionXS(float fission_xs) {
	_fission_xs = new float[1];
	_fission_xs_energies = new float[1];
	_fission_xs[0] = fission_xs;
	_fission_xs_energies[0] = 1.0;
	_num_fission_xs = 1;
    float (Isotope::*func)(float) const;
    func = &Isotope::getFissionXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    											const>(FISSION, func));
}


/**
 * Set the total cross-section for this isotope
 * @param total_xs a float array of microscopic total xs (barns)
 * @param total_xs_energies a float array of energies (eV)
 * @param num_total_xs the number of total xs values
 */
void Isotope::setTotalXS(float* total_xs, float* total_xs_energies,
													int num_total_xs) {

	_total_xs = total_xs;
    _total_xs_energies = total_xs_energies;
    _num_total_xs = num_total_xs;
    float (Isotope::*func)(float) const;
    func = &Isotope::getTotalXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
   													const>(TOTAL, func));
}


/**
 * Set the one total elastic scattering cross-section for this isotope
 * @param total_xs the one group total cross-section (barns)
 */
void Isotope::setOneGroupTotalXS(float total_xs) {
	_total_xs = new float[1];
	_total_xs_energies = new float[1];
	_total_xs[0] = total_xs;
	_total_xs_energies[0] = 1.0;
	_num_total_xs = 1;
    float (Isotope::*func)(float) const;
    func = &Isotope::getTotalXS;
    _xs_handles.insert(std::pair<collisionType, float(Isotope::*)(float)
    												const>(TOTAL, func));
}


/**
 * Add an inelastic scattering cdf to this isotope which applies to incoming
 * neutron energies between a lower and upper bound
 * @param lower_bound the lower bound incoming neutron energy (eV)
 * @param upper_bound the upper bound incoming neutron energy (eV)
 * @param energies a float array of outgoing energies
 * @param prob a float array of probabilities
 * @param num_bins the number of energy bins
 */
void Isotope::addInelasticScatterCDF(float lower_bound, float upper_bound,
					float* energies, float* prob, int num_bins) {

	/* Create a pair of the energy and probability arrays */
	std::pair<float*, float*> new_cdf = std::pair<float*, float*>
													(energies, prob);

	/* Create a pair of the previous pair along with the number of bins */
	std::pair<std::pair<float*, float*>, int> new_cdf_binned =
					std::pair<std::pair<float*, float*>, int>(new_cdf,
																num_bins);

	/* Create a pair of the upper and lower bound incoming neutron energies */
	std::pair<float, float> new_cdf_bounds =
					std::pair<float, float>(lower_bound, upper_bound);

	/* Insert this mess of pairs into the inelastic scattering cdf map */
	_inelastic_cdfs.insert(std::pair<std::pair<float, float>,
			std::pair<std::pair<float*, float*>, int> >(new_cdf_bounds,
														new_cdf_binned));
}


/**
 * For a given energy, this method determines a random collision type
 * based on the values of each of its cross-section types at that energy
 * @param energy the incoming neutron energy (eV)
 * @return the collision type (CAPTURE, SCATTER, INELASTIC, ELASTIC, TOTAL)
 */
collisionType Isotope::getCollisionType(float energy) {

	float test = float(rand()) / RAND_MAX;
	float collision_xs = 0.0;
	float next_collision_xs = 0.0;
	float total_xs = getTotalXS(energy);
	collisionType type = TOTAL;

	/* Loops over all cross-section types to find the one for this energy */
	std::map<collisionType, float(Isotope::*)(float)
												const>::const_iterator iter;

	for (iter = _xs_handles.begin(); iter != _xs_handles.end(); ++iter) {
		next_collision_xs += (this->*iter->second)(energy) / total_xs;

		if (test >= collision_xs && test <= next_collision_xs) {
			type = iter->first;
			break;
		}

		/* Update the next collision xs */
		collision_xs = next_collision_xs;
	}

	return type;
}


/**
 * For a given energy, this method determines a random collision type
 * based on the values of each of its cross-section types at that energy
 * @param energy the incoming neutron energy (eV)
 * @return the collision type (CAPTURE, SCATTER, INELASTIC, ELASTIC, TOTAL)
 */
collisionType Isotope::getOneGroupCollisionType() {

	float test = float(rand()) / RAND_MAX;
	float collision_xs = 0.0;
	float total_xs = getOneGroupTotalXS();
	collisionType type = TOTAL;

	/* Loops over all cross-section types to find the one for this energy */
	std::map<collisionType, float(Isotope::*)(float)
												const>::const_iterator iter;

	for (iter = _xs_handles.begin(); iter != _xs_handles.end(); ++iter) {
		collision_xs += (this->*iter->second)(1.0) / total_xs;

		if (test <= collision_xs)
			type = iter->first;
	}

	return type;
}


/**
 * For a given neutron energy in eV in a scattering collision, this
 * function returns the outgoing energy in eV, Eprime, for the collision
 * based on the thermal scattering distributions
 * @param energy the incoming energy (eV)
 * @return the outgoing energy (eV)
 */
float Isotope::getThermalScatteringEnergy(float energy) {

	/* First check that the thermal scattering CDFs have been initialized */
	if (_num_thermal_cdfs == 0)
		log_printf(ERROR, "Unable to sample the thermal scattering CDFs for"
				" isotope %s because they have not yet been initialized",
																_isotope_name);

	/* Convert energies in eV to eV / kT */
	energy /= (_kB * _T);

	/* Compute possible values for E to scatter to */
	float* possible_Eprimes = new float[_num_thermal_cdf_bins];
	for (int i=0; i < _num_thermal_cdf_bins; i++)
		possible_Eprimes[i] = _Eprime_to_E[i] * energy;

	float rn = float(rand()) / RAND_MAX;
	int index;
	float Eprime;

	/* Check if energy is lower than all thermal scattering CDFs */
	if (energy < _E_to_kT[0]) {
		index = findUpperIndex(_thermal_cdfs[0],
								_num_thermal_cdf_bins-1, 0, rn);
		Eprime = possible_Eprimes[index];
	}

	/* Check if energy is above all thermal scattering CDFs */
	else if (energy > _E_to_kT[_num_thermal_cdfs-1]) {
		index = findUpperIndex(_thermal_cdfs[_num_thermal_cdfs-1],
									_num_thermal_cdf_bins-1, 0, rn);
		Eprime = possible_Eprimes[index];
	}

	/* Otherwise the energy is sandwiched within the scattering CDFs */
	else {
		int upper_index = findUpperIndex(_E_to_kT, _num_thermal_cdfs-1,
															0, energy);
		int lower_index = upper_index - 1;
		int Eprime_lower_index = findUpperIndex(_thermal_cdfs[lower_index],
										_num_thermal_cdf_bins-1, 0, rn);
		float Eprime_lower = possible_Eprimes[Eprime_lower_index];
		int Eprime_upper_index = findUpperIndex(_thermal_cdfs[upper_index],
										_num_thermal_cdf_bins-1, 0, rn);
		float Eprime_upper = possible_Eprimes[Eprime_upper_index];
		float delta_E_to_kT = _E_to_kT[upper_index] - _E_to_kT[lower_index];
		float delta_Eprime = Eprime_upper - Eprime_lower;
		float slope = delta_Eprime / delta_E_to_kT;
		Eprime = slope * (energy - _E_to_kT[lower_index]) + Eprime_lower;
	}

	/* Convert outgoing energy back into eV */
	Eprime *= (_kB * _T);

	delete [] possible_Eprimes;

	return Eprime;
}


/**
 * This method clones a given Isotope class object by executing a deep
 * copy of all of the Isotope's class attributes and giving them to a new
 * Isotope class object
 * @return a pointer to the new cloned Isotope class object
 */
Isotope* Isotope::clone() {

	/* Allocate memory for the clone */
	Isotope* new_clone = new Isotope();

	/* Set the clones isotope name, atomic number, number density */
	new_clone->setIsotopeType(_isotope_name);
	new_clone->setA(_A);
	new_clone->setN(_N);
	new_clone->setTemperature(_T);

	/* If the given isotope has a one group capture xs */
	if (_num_capture_xs == 1)
		new_clone->setOneGroupCaptureXS(_capture_xs[0]);


	/* If the given isotope has an capture xs */
	if (_num_capture_xs > 1) {

		/* Deep copy the xs values */
		float* capture_xs = new float[_num_capture_xs];
		memcpy(capture_xs, _capture_xs, sizeof(float)*_num_capture_xs);

		/* Deep copy the energies for each of the xs values */
		float* capture_xs_energies = new float[_num_capture_xs];
		memcpy(capture_xs_energies, _capture_xs_energies,
				sizeof(float)*_num_capture_xs);

		/* Set the clone's capture xs */
		new_clone->setCaptureXS(capture_xs, capture_xs_energies, _num_capture_xs);
	}


	/* If the given isotope has a one group scatter xs */
	if (_num_scatter_xs == 1)
		new_clone->setOneGroupScatterXS(_scatter_xs[0], _scatter_angle);


	/* If the given isotope has a scatter xs */
	if (_num_scatter_xs > 1) {

		/* Deep copy the xs values */
		float* scatter_xs = new float[_num_scatter_xs];
		memcpy(scatter_xs, _scatter_xs, sizeof(float)*_num_scatter_xs);

		/* Deep copy the energies for each of the xs values */
		float* scatter_xs_energies = new float[_num_scatter_xs];
		memcpy(scatter_xs_energies, _scatter_xs_energies,
								sizeof(float)*_num_scatter_xs);

		/* Set the clone's xs */
		new_clone->setScatterXS(scatter_xs, scatter_xs_energies,
				_num_scatter_xs, _scatter_angle);
	}


	/* If the given isotope has an inelastic scatter xs */
	if (_num_inelastic_xs == 1)
		new_clone->setOneGroupInelasticXS(_inelastic_xs[0], _inelastic_angle);


	/* If the given isotope has an inelastic scatter xs */
	if (_num_inelastic_xs > 1) {

		/* Deep copy the xs values */
		float* inelastic_xs = new float[_num_inelastic_xs];
		memcpy(inelastic_xs, _inelastic_xs,
				sizeof(float)*_num_inelastic_xs);

		/* Deep copy the energies for each of the xs values */
		float* inelastic_xs_energies = new float[_num_inelastic_xs];
		memcpy(inelastic_xs_energies, _inelastic_xs_energies,
							sizeof(float)*_num_inelastic_xs);

		/* Set the clone's xs */
		new_clone->setInelasticXS(inelastic_xs, inelastic_xs_energies,
				_num_inelastic_xs, _inelastic_angle);
	}


	/* If the given isotope has an elastic scatter xs */
	if (_num_elastic_xs == 1)
		new_clone->setOneGroupElasticXS(_elastic_xs[0], _elastic_angle);


	/* If the given isotope has an elastic scatter xs */
	if (_num_elastic_xs > 1) {

		/* Deep copy the xs values */
		float* elastic_xs = new float[_num_elastic_xs];
		memcpy(elastic_xs, _elastic_xs, sizeof(float)*_num_elastic_xs);

		/* Deep copy the energies for each of the xs values */
		float* elastic_xs_energies = new float[_num_elastic_xs];
		memcpy(elastic_xs_energies, _elastic_xs_energies,
								sizeof(float)*_num_elastic_xs);

		/* Set the clone's xs */
		new_clone->setElasticXS(elastic_xs, elastic_xs_energies,
				_num_elastic_xs, _elastic_angle);
	}

	/* If the given isotope has a one group fission xs */
	if (_num_fission_xs == 1)
		new_clone->setOneGroupFissionXS(_fission_xs[0]);


	/* If the given isotope has a fission xs */
	if (_num_fission_xs > 1) {

		/* Deep copy the xs values */
		float* fission_xs = new float[_num_fission_xs];
		memcpy(fission_xs, _fission_xs, sizeof(float)*_num_fission_xs);

		/* Deep copy the energies for each of the xs values */
		float* fission_xs_energies = new float[_num_fission_xs];
		memcpy(fission_xs_energies, _fission_xs_energies,
				sizeof(float)*_num_fission_xs);

		/* Set the clone's fission xs */
		new_clone->setFissionXS(fission_xs, fission_xs_energies,
													_num_fission_xs);
	}


	/* If the given isotope has an explicitly defined total scatter xs */
	if (_num_total_xs == 1)
		new_clone->setOneGroupTotalXS(_total_xs[0]);


	/* If the given isotope has an explicitly defined total scatter xs */
	if (_num_total_xs > 1) {

		/* Deep copy the xs values */
		float* total_xs = new float[_num_total_xs];
		memcpy(total_xs, _total_xs, sizeof(float)*_num_total_xs);

		/* Deep copy the energies for each of the xs values */
		float* total_xs_energies = new float[_num_total_xs];
		memcpy(total_xs_energies, _total_xs_energies,
						sizeof(float)*_num_total_xs);

		/* Set the clone's xs */
		new_clone->setTotalXS(total_xs, total_xs_energies,
				_num_total_xs);
	}

	/* Deep copy the isotope's inelastic scatter cdfs */
	std::map<std::pair<float, float>, std::pair< std::pair<float*,
									float*>, int> >::iterator iter;

	for (iter = _inelastic_cdfs.begin(); iter!=_inelastic_cdfs.end(); ++iter) {

		float* cdf_energies = new float[iter->second.second];
		memcpy(cdf_energies, iter->second.first.first,
							sizeof(float)*iter->second.second);
		float* prob = new float[iter->second.second];
		memcpy(prob, iter->second.first.second,
							sizeof(float)*iter->second.second);

		new_clone->addInelasticScatterCDF(iter->first.first,
				iter->first.second, cdf_energies, prob, iter->second.second);

	}

	/* Initialize the isotope's thermal scattering CDFs if they have been
	 * created for this isotope */
	if (_num_thermal_cdfs > 0)
		new_clone->initializeThermalScattering(_E_to_kT[0]*_kB*_T,
			_E_to_kT[_num_thermal_cdfs-1]*_kB*_T, _num_thermal_cdf_bins,
													_num_thermal_cdfs);

	/* Return a pointer to the cloned Isotope class */
	return new_clone;
}



/* This method initializes the probability distributions for thermal
 * scattering. It takes in arguments for the start_energy and end
 * energy (ratios of kT) and the number of distributions which it
 * uses to generate logarithmically spaced energies for the distributions.
 * It also takes in the number of energy bins to use for each bin.
 * @param start_energy the first distribution's energy
 * @param end_energy the final distribution's energy
 * @param num_bins the number of bins per distribution
 * @param end_distributions the number of scattering distributions
 */
void Isotope::initializeThermalScattering(float start_energy,
					float end_energy, int num_bins, int num_distributions) {

	/* Number of thermal scattering distributions */
	_num_thermal_cdfs = num_distributions;

	/* Number of bins per distribution */
	_num_thermal_cdf_bins = num_bins;

	/* Allocate memory for distributions */
	_thermal_cdfs = new float*[_num_thermal_cdfs];
	for (int i=0; i < _num_thermal_cdfs; i++)
		_thermal_cdfs[i] = new float[_num_thermal_cdf_bins];

	_thermal_dist = new float[_num_thermal_cdfs * _num_thermal_cdf_bins];
	float* cdf = new float[_num_thermal_cdf_bins];

	/* Initialize logarithmically spaced E/kT for each distribution */
	_E_to_kT = logspace(start_energy/(_kB*_T), end_energy/(_kB*_T),
												_num_thermal_cdfs);

	/* Find the maximum Eprime / E value that we must extend our distributions
	 * to before they all fall below some tolerance */
	bool tolerance_met = false;
	float dist_tolerance = 0.1;
	float curr_prob;
	float curr_Eprime_to_E = 1.0;
	while (!tolerance_met) {

		/* Start with the tolerance being met */
		tolerance_met = true;

		/* Loop over all CDFs and check if we are within the threshold */
		for (int i=0; i < _num_thermal_cdfs; i++) {
			curr_prob = thermalScatteringProb(curr_Eprime_to_E, i);

			/* If we are above the tolerance */
			if (curr_prob > dist_tolerance)
				tolerance_met = false;
		}

		/* Update distance along x-axis */
		if(!tolerance_met)
			curr_Eprime_to_E += 0.25;
	}

	/* Initialize x-axis of Eprime to E ratios */
	_Eprime_to_E = logspace(1E-5, curr_Eprime_to_E, _num_thermal_cdf_bins+1);

	/* Loop over each distribution */
	for (int i=0; i < _num_thermal_cdfs; i++) {
		for (int j=0; j < _num_thermal_cdf_bins; j++)
			_thermal_dist[i*_num_thermal_cdf_bins + j] =
									thermalScatteringProb(_Eprime_to_E[j], i);
	}

	/* Create CDFs for each distribution */
	for (int i=0; i < _num_thermal_cdfs; i++) {
		cumulativeIntegral(_Eprime_to_E,
							&_thermal_dist[i*_num_thermal_cdf_bins], cdf,
										_num_thermal_cdf_bins, TRAPEZOIDAL);

		/* Transfer CDF values to our array */
		for (int j=0; j < _num_thermal_cdf_bins; j++)
			_thermal_cdfs[i][j] = cdf[j];
	}    void rescaleXS(float* new_energies, int num_energies);


	delete [] cdf;

	/* Normalize CDFs */
	for (int i=0; i < _num_thermal_cdfs; i++) {
		for (int j=0; j < _num_thermal_cdf_bins; j++)
			_thermal_cdfs[i][j] /= _thermal_cdfs[i][_num_thermal_cdf_bins-1];
	}

	return;
}


/**
 * This function computes the thermal scattering probability for
 * for a ratio of initial to final energies
 * @param E_prime_to_E a ratio of initial to final energies
 * @return the probability of the ratio occurring
 */
float Isotope::thermalScatteringProb(float E_prime_to_E, int dist_index) {

	float prob;

    /* Computes the final energy for each of the ratios */
    float Eprime = _E_to_kT[dist_index] * E_prime_to_E;

    /* Uses the equation from 22.211 slide 26 of the 2nd lecture
     * to compute probabilities */
    float a = sqrt(_E_to_kT[dist_index]);
	float b = sqrt(Eprime);
	float c = erf(_eta * b - _rho * a);
	float d = erf(_eta * b + _rho * a);
	float e = erf(_eta * a - _rho * b);
	float f = erf(_eta * a + _rho * b);
	double g = exp(double(_E_to_kT[dist_index]) - double(Eprime));

	/* Account for lower and upper signs in equation */
	if (Eprime > _E_to_kT[dist_index])
		prob = (c - d) + g * (e + f);
	else
		prob = (c + d) + g * (e - f);

	/* Multiply by eta / 2 */
	prob *= _eta*_eta / 2.0;

	/* Normalize to the atomic mass by multiplying by 1-alpha */
	prob *= (1.0 - _alpha);

	return prob;
}



void Isotope::rescaleXS(float* energies, int num_energies) {

	/* Loops over all cross-section types to find the one for this energy */
	std::map<collisionType, float(Isotope::*)(float)
												const>::const_iterator iter;

	for (iter = _xs_handles.begin(); iter != _xs_handles.end(); ++iter) {

		float* new_xs = new float[num_energies];
		float* new_energies = new float[num_energies];
		memcpy(new_energies, energies, sizeof(float)*num_energies);

		for (int i=0; i < num_energies; i++)
			new_xs[i] = (this->*iter->second)(new_energies[i]);

		if (iter->first == CAPTURE) {
			_num_capture_xs = num_energies;
			delete [] _capture_xs_energies;
			delete _capture_xs;
			setCaptureXS(new_xs, new_energies, num_energies);
		}
		else if (iter->first == SCATTER) {
			_num_scatter_xs = num_energies;
			delete [] _scatter_xs_energies;
			delete _scatter_xs;
			setScatterXS(new_xs, new_energies, num_energies, _scatter_angle);
		}
		else if (iter->first == ELASTIC) {
			_num_elastic_xs = num_energies;
			delete [] _elastic_xs_energies;
			delete _elastic_xs;
			setElasticXS(new_xs, new_energies, num_energies, _elastic_angle);
		}
		else if (iter->first == INELASTIC) {
			_num_inelastic_xs = num_energies;
			delete [] _inelastic_xs_energies;
			delete _inelastic_xs;
			setInelasticXS(new_xs, new_energies, num_energies, _inelastic_angle);
		}
		else if (iter->first == FISSION) {
			_num_fission_xs = num_energies;
			delete [] _fission_xs_energies;
			delete _fission_xs;
			setFissionXS(new_xs, new_energies, num_energies);
		}
		else if (iter->first == TOTAL) {
			_num_total_xs = num_energies;
			delete [] _total_xs_energies;
			delete _total_xs;
			setTotalXS(new_xs, new_energies, num_energies);
		}
	}

	return;
}

/**
 * This method plots the microscopic cross-section values for variable number
 * of cross-section types for this isotope. The cross-section types are
 * specified followed by the number "-1" after all types have been listed
 * isotope
 * @param start_energy the energy to begin the plot at
 * @param end_energy the energy to end the plot at
 * @param num_energies the number of energies to plot
 * @param type the cross-section types to plot
 */
void Isotope::plotXS(float start_energy, float end_energy, int num_energies,
													collisionType types, ...) {

	/* Create an array of logarithmically spaced energies */
	float* energies = logspace(start_energy, end_energy, num_energies);
	float* xs_values = new float[num_energies];

	/* Initialize variable parameters data structures of different isotopes */
	va_list xs_types;
	va_start(xs_types, types);
	int i;
    float (Isotope::*func)(float) const;
    func = &Isotope::getTotalXS;

	/* Create title and filename for plot */
	std::stringstream title;
	std::stringstream filename;
	title << "set title \"" << _isotope_name;
	title << " Microscopic Cross-sections\"";
	filename << _isotope_name << "_micro_xs";

	/* Initialize the plot */
	gnuplot_ctrl* handle = gnuplot_init();
	gnuplot_set_xlabel(handle, (char*)"Energy (eV)");
	gnuplot_set_ylabel(handle, (char*)"Cross-section (barns)");
	gnuplot_cmd(handle, (char*)title.str().c_str());
	gnuplot_cmd(handle, (char*)"set logscale xy");
	gnuplot_setstyle(handle, (char*)"lines");

	/* Loop through each isotope */
	for (i=types; i >= 0; i=va_arg(xs_types, int)) {

		std::stringstream legend;

		/* Get the appropriate function handle to compute the cross-section
		 * type the user requested */
		if (i == CAPTURE) {
			legend << "capture";
			func = &Isotope::getCaptureXS;
		}
		else if (i == ELASTIC) {
			legend << "elastic";
			func = &Isotope::getElasticXS;
		}
		else if (i == INELASTIC) {
			legend << "inelastic";
			func = &Isotope::getInelasticXS;
		}
		else if (i == SCATTER) {
			legend << "scatter";
			func = &Isotope::getScatterXS;
		}
		else if (i == FISSION) {
			legend << "fission";
			func = &Isotope::getFissionXS;
		}

	    /* Compute the cross-section at each of the energies */
		for (int i=0; i < num_energies; i++)
			xs_values[i] = (this->*func)(energies[i]);

		/* Plot the cross-section */
		gnuplot_saveplot(handle, (char*)filename.str().c_str());
		gnuplot_plot_xy(handle, energies, xs_values, num_energies,
									(char*)legend.str().c_str());
	}

	gnuplot_close(handle);
	va_end(xs_types);

	delete [] energies;
	delete [] xs_values;
}


/**
 * Plot the thermal scattering distributions for this isotope if they have
 * been initialized
 */
void Isotope::plotThermalScatteringDistributions() {

	/* Check that the thermal scattering cdfs have been initialized */
	if (_num_thermal_cdfs == 0)
		log_printf(ERROR, "Unable to plot the thermal scattering distributions"
				" because they have not yet been initialized");

	std::stringstream title;
	std::stringstream filename;
	title << "set title \"" << _isotope_name << " Thermal Scattering Dist\"";
	filename << _isotope_name << "_therm_dist";

	/* Plot the Cross-section */
	gnuplot_ctrl* handle = gnuplot_init();
	gnuplot_set_xlabel(handle, (char*)"Eprime / E");
	gnuplot_set_ylabel(handle, (char*)"PDF");
	gnuplot_cmd(handle, (char*)title.str().c_str());
	gnuplot_setstyle(handle, (char*)"lines");
	for (int i=0; i < _num_thermal_cdfs; i++) {
		std::stringstream legend;
		legend << _E_to_kT[i] << "kT";
		if (i == _num_thermal_cdfs - 1)
			gnuplot_saveplot(handle, (char*)filename.str().c_str());
		gnuplot_plot_xy(handle, _Eprime_to_E,
					&_thermal_dist[i*_num_thermal_cdf_bins],
					_num_thermal_cdf_bins-1, (char*)legend.str().c_str());
	}
	gnuplot_close(handle);

}


/**
 * Plot the thermal scattering CDFs for this isotope type
 */
void Isotope::plotThermalScatteringCDFs() {

	/* Check that the thermal scattering cdfs have been initialized */
	if (_num_thermal_cdfs == 0)
		log_printf(ERROR, "Unable to plot the thermal scattering CDFs because"
				"they have not yet been initialized");

	std::stringstream title;
	std::stringstream filename;
	title << "set title \"" << _isotope_name << " Thermal Scattering CDFs\"";
	filename << _isotope_name << "_therm_cdfs";

	/* Plot the Cross-section */
	gnuplot_ctrl* handle = gnuplot_init();
	gnuplot_set_xlabel(handle, (char*)"Eprime / E");
	gnuplot_set_ylabel(handle, (char*)"CDF");
	gnuplot_cmd(handle, (char*)title.str().c_str());
	gnuplot_setstyle(handle, (char*)"lines");
	for (int i=0; i < _num_thermal_cdfs; i++) {
		std::stringstream legend;
		legend << _E_to_kT[i] << "kT";
		if (i == _num_thermal_cdfs - 1)
			gnuplot_saveplot(handle, (char*)filename.str().c_str());
		gnuplot_plot_xy(handle, _Eprime_to_E, _thermal_cdfs[i],
					_num_thermal_cdf_bins-1, (char*)legend.str().c_str());
	}
	gnuplot_close(handle);
}


/**
 * This method plots the isotope's thermal scattering CDFs sampled at a
 * certain energy
 * @param energy the energy (eV) to sample
 * @param num_samples the number of samples to plot
 */
void Isotope::plotSampledThermalScatteringEnergies(float energy,
												int num_samples) {

	/* Create a new Binner to bin the sampled outgoing energies */
	Binner* bins = new Binner();
	bins->generateBinEdges(0.0, energy*30, 50, EQUAL);

	float sample;

	/* Generate sampled outgoing energies */
	for (int i=0; i < num_samples; i++) {
		sample = getThermalScatteringEnergy(energy);
		bins->tally(sample);
	}

	std::stringstream title;
	std::stringstream filename;
	title << "set title \"" << _isotope_name << " Sampled Thermal Scattering "
			<< "Eprimes (for E = " << energy << ")\"";
	filename << _isotope_name << "_" << energy << "eV_therm_cdfs";

	/* Plot the Cross-section */
	gnuplot_ctrl* handle = gnuplot_init();
	gnuplot_set_xlabel(handle, (char*)"Eprime (eV)");
	gnuplot_set_ylabel(handle, (char*)"PDF");
	gnuplot_cmd(handle, (char*)title.str().c_str());
	gnuplot_setstyle(handle, (char*)"lines");
	gnuplot_saveplot(handle, (char*)filename.str().c_str());
	gnuplot_plot_xy(handle, bins->getBinCenters(), bins->getTallies(),
								bins->getNumBins(), (char*)"Samples");
}
