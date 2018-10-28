#ifndef _MYWORLD_
#define _MYWORLD_

#include <vector>
#include <iostream>

#define IX(i, j) ((i)+(getNumCells()+2)*(j))
#define IXColor(i, j, c) ((i)*(getNumChannels())*(getNumCells()+2)+(j)*(getNumChannels())+(c))
#define SWAPPING(x0,x) {double *tmp=x0;x0=x;x=tmp;}
#define IXRed(i, j) IXColor(i, j, 0)
#define IXGreen(i,j) IXColor(i, j, 1)
#define IXBlue(i,j) IXColor(i, j, 2)

class MyWorld {
 public:
    MyWorld();

    virtual ~MyWorld();

    void initialize(int _numCells, double _timeStep, double _diffCoef, double _viscCoef, int channels);
    double getTimeStep();
    int getNumCells() { return mNumCells;}
    double getDensity(int _index) { return mDensity[_index]; }
    double getVelocityU(int _index) { return mU[_index]; }
    double getVelocityV(int _index) { return mV[_index]; }
    void setDensity(int _i, int _j, double _source) { mDensity[IX(_i, _j)] += mTimeStep * _source; }
    void setDensity(int _i, int _j, int c, double _source) { mDensity[IXColor(_i, _j, c)] += mTimeStep * _source; }
    void setDensity(int _i, double _source) { mDensity[_i] = _source;}
    void setU(int _i, int _j, double _force) { mU[IX(_i, _j)] += mTimeStep * _force; }
    void setV(int _i, int _j, double _force) { mV[IX(_i, _j)] += mTimeStep * _force; }
    void reset();
    void simulate();
    int getNumChannels() { return mChannels;}
    void printDensity() { for(int i = 0; i < (2*getNumCells()*getNumChannels()); i++){ std::cout << std::hex << (int)mDensity[i] <<"," << std::endl; }}
    unsigned char* getDensityAsByte();
 protected:
    void densityStep(double *_x, double *_x0);
    void velocityStep(double *_u, double *_v, double *_u0, double *_v0);
    void diffuseDensity(double *_x, double *_x0);
    void diffuseVelocity(double *_u, double *_v, double *_u0, double *_v0);
    void advect(double *_d, double *_d0, double *_u, double *_v, int channels=1);
    void advectDensity(double *_d, double *_d0, double *_u, double *_v);
    void advectVelocity(double *_u, double *_v, double *_u0, double *_v0);
    void project(double *_u, double *_v, double *_u0, double *_v0);
    void externalForces();
    void linearSolve(double *_x, double *_x0, double _a, double _c, int channels=1);
    void setBoundary(double *_x);
    void setVelocityBoundary(double *_u, double *_v);
    void setBoundaryWithChannels(double *_x);
    void setBoundaryChannelsHelper(double *_x, int c);


    int mChannels;
    int mNumCells;
    double mTimeStep;
    double mDiffusionCoef;
    double mViscosityCoef;
    double *mU;
    double *mV;
    double *mPreU;
    double *mPreV;
    double *mDensity;
    unsigned char *mDensityInBytes;
    double *mPreDensity;
};

#endif
