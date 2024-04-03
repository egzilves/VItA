/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CCOCommonStructures.h
 *
 *  Created on: Jun 1, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef CCOCOMMONSTRUCTURES_H_
#define CCOCOMMONSTRUCTURES_H_

#include <vector>

#include <vtkSmartPointer.h>
#include <vtkLine.h>

using namespace std;

typedef double (*nuBloodFunc)(int treeLevel);
typedef double (*gammaFunc)(int treeLevel);
typedef double (*epsLimFunc)(int treeLevel);

/**
 * Point structure endowed with basic linear algebra operations.
 */
struct point {
	double p[3];

	/**	Element to element addition.*/
	inline point operator+(point a) {
		return {a.p[0]+p[0],a.p[1]+p[1], a.p[2]+p[2]};
	}
	/**	Element to element product.*/
	inline point operator*(point a) {
		return {a.p[0]*p[0],a.p[1]*p[1], a.p[2]*p[2]};
	}
	/**	Scalar to element division.*/
	inline point operator/(double alpha) {
		return {p[0]/alpha,p[1]/alpha, p[2]/alpha};
	}
	/**	Scalar to element product.*/
	inline point operator*(double alpha) {
		return {alpha*p[0],alpha*p[1], alpha*p[2]};
	}
	/**	Inner product */
	inline double operator^(point a) {
		return a.p[0] * p[0] + a.p[1] * p[1] + a.p[2] * p[2];
	}
	/**	Element to element substraction.*/
	inline point operator-(point a) {
		return {p[0]-a.p[0],p[1]-a.p[1], p[2]-a.p[2]};
	}
	/**	Element to element equallity.*/
	inline int operator==(point a) {
		return (a.p[0] == p[0] && a.p[1] == p[1] && a.p[2] == p[2]);
	}
	// /**	Outter product */
	// inline struct matrix operator%(point a) {
	// 	return {
	// 			p[0]*a.p[0], p[0]*a.p[1], p[0]*a.p[2], 
	// 			p[1]*a.p[0], p[1]*a.p[1], p[1]*a.p[2], 
	// 			p[2]*a.p[0], p[2]*a.p[1], p[2]*a.p[2]
	// 		};
	// }
};
/**
 * Matrix structure endowed with basic linear algebra operations.
 */
struct matrix {
	double e[3][3];

	/**	Element to element addition.*/
	inline matrix operator+(matrix a) {
		return {
			e[0][0] + a.e[0][0], e[0][1] + a.e[0][1], e[0][2] + a.e[0][2],
			e[1][0] + a.e[1][0], e[1][1] + a.e[1][1], e[1][2] + a.e[1][2],
			e[2][0] + a.e[2][0], e[2][1] + a.e[2][1], e[2][2] + a.e[2][2]
		};
	}
	/**	Element to element subtraction.*/
	inline matrix operator-(matrix a) {
		return {
			e[0][0] - a.e[0][0], e[0][1] - a.e[0][1], e[0][2] - a.e[0][2],
			e[1][0] - a.e[1][0], e[1][1] - a.e[1][1], e[1][2] - a.e[1][2],
			e[2][0] - a.e[2][0], e[2][1] - a.e[2][1], e[2][2] - a.e[2][2]
		};
	}
	/**	Scalar to element product.*/
	inline matrix operator*(double alpha) {
		return {
			alpha*e[0][0], alpha*e[0][1], alpha*e[0][2],
			alpha*e[1][0], alpha*e[1][1], alpha*e[1][2],
			alpha*e[2][0], alpha*e[2][1], alpha*e[2][2]
		};
	}
	/**	Scalar to element division.*/
	inline matrix operator/(double alpha) {
		return {
			e[0][0]/alpha, e[0][1]/alpha, e[0][2]/alpha,
			e[1][0]/alpha, e[1][1]/alpha, e[1][2]/alpha,
			e[2][0]/alpha, e[2][1]/alpha, e[2][2]/alpha
		};
	}
	/** Matrix multiplication.*/
	inline matrix operator*(matrix a){
		return{
			e[0][0]*a.e[0][0] + e[0][1]*a.e[1][0] + e[0][2]*a.e[2][0], e[0][0]*a.e[0][1] + e[0][1]*a.e[1][1] + e[0][2]*a.e[2][1], e[0][0]*a.e[0][2] + e[0][1]*a.e[1][2] + e[0][2]*a.e[2][2], 
			e[1][0]*a.e[0][0] + e[1][1]*a.e[1][0] + e[1][2]*a.e[2][0], e[1][0]*a.e[0][1] + e[1][1]*a.e[1][1] + e[1][2]*a.e[2][1], e[1][0]*a.e[0][2] + e[1][1]*a.e[1][2] + e[1][2]*a.e[2][2], 
			e[2][0]*a.e[0][0] + e[2][1]*a.e[1][0] + e[2][2]*a.e[2][0], e[2][0]*a.e[0][1] + e[2][1]*a.e[1][1] + e[2][2]*a.e[2][1], e[2][0]*a.e[0][2] + e[2][1]*a.e[1][2] + e[2][2]*a.e[2][2]
		};
	}
	/** Matrix-vector multiplication.*/
	inline matrix operator*(point a){
		return{
			e[0][0]*a.p[0] + e[0][1]*a.p[1] + e[0][2]*a.p[2],
			e[1][0]*a.p[0] + e[1][1]*a.p[1] + e[1][2]*a.p[2],
			e[2][0]*a.p[0] + e[2][1]*a.p[1] + e[2][2]*a.p[2]
		};
	}
	/**	Inner product */
	inline double operator^(matrix a) {
		return {a.e[0][0] * e[0][0] + a.e[0][1] * e[0][1] + a.e[0][2] * e[0][2] + 
			a.e[1][0] * e[1][0] + a.e[1][1] * e[1][1] + a.e[1][2] * e[1][2] + 
			a.e[2][0] * e[2][0] + a.e[2][1] * e[2][1] + a.e[2][2] * e[2][2]};
	}
};
inline matrix outer(point a, point b){
	return matrix {
		a.p[0]*b.p[0], a.p[0]*b.p[1], a.p[0]*b.p[2], 
		a.p[1]*b.p[0], a.p[1]*b.p[1], a.p[1]*b.p[2], 
		a.p[2]*b.p[0], a.p[2]*b.p[1], a.p[2]*b.p[2]
	};
}

/**
 * Overload operator to print point @p p in a stream object.
 * @param os	Stream output.
 * @param p		Point to print.
 * @return	Streamed output.
 */
inline ostream& operator<<(ostream& os, point p) {
	os << "(" << p.p[0] << ", " << p.p[1] << ", " << p.p[2] << ")";
	return os;
}

/**
 * Vessel structure for CCO trees.
 */
struct vessel {
	/**
	 * Anastomose of the vessel. Position 0 is the parent vessel.
	 */
	vector<vessel *> anastomose;

	/**
	 * VTK geomtric representation.
	 */
	vtkSmartPointer<vtkLine> vtkSegment;
	/**
	 * Unique identifier for the vessel.
	 */
	vtkIdType vtkSegmentId;

	/**
	 * Proximal position of the vessel.
	 */
	point xProx;
	/**
	 * Distal position of the vessel.
	 */
	point xDist;

	/** Bifurcation level from proximal to distal (Root is at level 0).*/
	int nLevel;
	/** Vessel radius. */
	double radius;
	/** Radius relative to the parent vessel (i.e. radius/parent_radius). For root, it is the radius value. */
	double beta;
	/** Distance between xProx and xDist. */
	double length;
	/** Reduced fluid-dynamic resistance. */
	double resistance;
	/** Reduced fluid-dynamic resistance. */
	double localResistance;
	/** Blood viscosity. */
	double viscosity;
	/** Flow in this vessel. */
	double flux;
	/** Pressure in the vessel. */
	double pressure;
	/** Generation step */
	long long int ID;
	/** Volume of this down tree branch.*/
	double treeVolume;
};

/**
 * Overload operator to print vessel @p v in a stream object.
 * @param os	Stream output.
 * @param v		Vessel to print.
 * @return	Streamed output.
 */
inline ostream& operator<<(ostream& os, vessel *v) {
	os << "Vessel " << v->vtkSegmentId << " at level " << v->nLevel << " with ptr=" << (void *) v << endl;
	os << "Location " << v->xProx << " - " << v->xDist << " with a length of " << v->length << " cm" << endl;
	os << "Beta=" << v->beta << ", Rsub=" << v->resistance << ", flux=" << v->flux << " cm^3/s and subtree cost of " << v->treeVolume << endl;
	os << "Anastomose (Ids): ";
	for (vector<vessel *>::iterator it = v->anastomose.begin(); it < v->anastomose.end(); ++it)
		if (*it)
			os << (*it)->vtkSegmentId << " ";
		else
			os << "-1";
	os << endl;
	return os;
}

#endif /* CCOCOMMONSTRUCTURES_H_ */
