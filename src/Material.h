/*
 * Material.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef MATERIAL_H_
#define MATERIAL_H_

#include <map>
#include <math.h>
#include <stdarg.h>
#include "interpolate.h"
#include "integrate.h"
#include "gnuplot.h"
#include "xsreader.h"
#include "Binner.h"

/* Types of collisions */
typedef enum collisionTypes{
	CAPTURE,
	SCATTER,
	INELASTIC,
	ELASTIC,
	FISSION,
	TOTAL
} collisionType;

/* Types of angular scattering distributions */
typedef enum scatterAngleTypes {
	ISOTROPIC_CM,
	ISOTROPIC_LAB
} scatterAngleType;


/**
 * The material class represents an isotope and all of its material properties
 * which are relevant to neutronics
 */
class Material
{
private:
	char* _isotope;
	int _A;
	float _alpha;
	float _eta;
	float _rho;
	float _N;
	float _T;
	float _mu_avg;
	float* _capture_xs;
	float* _capture_xs_energies;
	int _num_capture_xs;
	float* _scatter_xs;
	float* _scatter_xs_energies;
	int _num_scatter_xs;
	scatterAngleType _scatter_angle;
	float* _inelastic_xs;
	float* _inelastic_xs_energies;
	int _num_inelastic_xs;
	scatterAngleType _inelastic_angle;
	float* _elastic_xs;
	float* _elastic_xs_energies;
	int _num_elastic_xs;
	scatterAngleType _elastic_angle;
	float* _fission_xs;
	float* _fission_xs_energies;
	int _num_fission_xs;
	float* _total_xs;
	float* _total_xs_energies;
	int _num_total_xs;

	/* Map of keys (xs types) with values (getXS functions for xs types) */
	std::map<collisionType, float(Material::*)(float) const> _xs_handles;

	/* Map of keys (incoming inelastic scatter cdf bounds) with values for
	 * inelastic scatter outgoing energy distributions */
	std::map<std::pair<float, float>, std::pair< std::pair<float*,
										float*>, int> > _inelastic_cdfs;

	int _num_thermal_cdfs;
	int _num_thermal_cdf_bins;
	float* _thermal_dist;
	float** _thermal_cdfs;
	float* _E_to_kT;
	float* _Eprime_to_E;
	float _kB;
public:
	Material();
    virtual ~Material();

    char* getIsotopeType() const;
    int getA() const;
    float getAlpha() const;
    float getN() const;
    float getTemperature() const;
    float getMuAverage() const;
    float getCaptureXS(float energy) const;
    float getOneGroupCaptureXS() const;
    float getScatterXS(float energy) const;
    float getOneGroupScatterXS() const;
    scatterAngleType getScatterAngleType() const;
    float getInelasticXS(float energy) const;
    float getOneGroupInelasticXS() const;
    scatterAngleType getInelasticAngleType() const;
    float getElasticXS(float energy) const;
    float getOneGroupElasticXS() const;
    scatterAngleType getElasticAngleType() const;
    float getFissionXS(float energy) const;
    float getOneGroupFissionXS() const;
    float getAbsorbXS(float energy) const;
    float getOneGroupAbsorbXS() const;
    float getTotalXS(float energy) const;
    float getOneGroupTotalXS() const;
    float getTransportXS(float energy) const;
    float getInelasticScatterEnergy(float energy);
    bool usesThermalScattering();

    void setIsotopeType(char* isotope);
    void setA(int A);
    void setN(float N);
    void setTemperature(float T);
    void loadXS(char* filename, collisionType type, char* delimiter);
    void setCaptureXS(float* capture_xs, float* capture_xs_energies,
    											int num_capture_xs);
    void setOneGroupCaptureXS(float capture_xs);
    void setScatterXS(float* scatter_xs, float* scatter_xs_energies,
							int num_scatter_xs, scatterAngleType type);
    void setOneGroupScatterXS(float scatter_xs, scatterAngleType type);
    void setScatterAngleType(scatterAngleType type);
    void setInelasticXS(float* inelastic_xs, float* inelastic_xs_energies,
							int num_inelastic_xs, scatterAngleType type);
    void setOneGroupInelasticXS(float inelastic_xs, scatterAngleType type);
    void setInelasticAngleType(scatterAngleType type);
    void setElasticXS(float* elastic_xs, float* elastic_xs_energies,
							int num_elastic_xs, scatterAngleType type);
    void setOneGroupElasticXS(float elastic_xs, scatterAngleType type);
    void setElasticAngleType(scatterAngleType type);
    void setFissionXS(float* fission_xs, float* fission_xs_energies,
												int num_fission_xs);
    void setOneGroupFissionXS(float fission_xs);
    void setTotalXS(float* scatter_xs, float* total_xs_energies,
												int num_total_xs);
    void setOneGroupTotalXS(float total_xs);
    void addInelasticScatterCDF(float lower_bound, float upper_bound,
    					float* energies, float* prob, int num_bins);

    collisionType getCollisionType(float energy);
    collisionType getOneGroupCollisionType();
    float getThermalScatteringEnergy(float energy);
    Material* clone();
    void initializeThermalScattering(float start_energy, float end_energy,
    								int num_bins, int num_distributions);
    float thermalScatteringProb(float E_prime_to_E, int dist_index);

    void plotXS(float start_energy, float end_energy, int num_energies,
    										collisionType type, ...);
    void plotThermalScatteringDistributions();
    void plotThermalScatteringCDFs();
    void plotSampledThermalScatteringEnergies(float energy, int num_samples);
};

#endif /* MATERIAL_H_ */
