#include "MyWorld.h"


MyWorld::MyWorld() {
}

MyWorld::~MyWorld() {
}

void MyWorld::initialize(int _numCells, double _timeStep, double _diffCoef, double _viscCoef, int channels) {
    mNumCells = _numCells;
    mTimeStep = _timeStep;
    mDiffusionCoef = _diffCoef;
    mViscosityCoef = _viscCoef;
    mChannels = channels;
    int size = (mNumCells + 2) * (mNumCells + 2);

    // Allocate memory for velocity and density fields
    mU = new double[size];
    mV = new double[size];
    mPreU = new double[size];
    mPreV = new double[size];
    mDensity = new double[size * mChannels];
    mPreDensity = new double[size * mChannels];

    reset();
}

void MyWorld::reset() {
    int size = (mNumCells + 2) * (mNumCells + 2);
    std::fill_n(mU, size, 0);
    std::fill_n(mV, size, 0);
    std::fill_n(mPreU, size, 0);
    std::fill_n(mPreV, size, 0);
    std::fill_n(mDensity, size * mChannels, 0);
    std::fill_n(mPreDensity, size * mChannels, 0);
}

double MyWorld::getTimeStep() {
    return mTimeStep;
}

void MyWorld::simulate() {
    velocityStep(mU, mV, mPreU, mPreV);
    densityStep(mDensity, mPreDensity);
    externalForces();
}

void MyWorld::densityStep(double *_x, double *_x0) {
    SWAPPING(_x, _x0); // _x now points at mPreDensity
    diffuseDensity(_x, _x0); // Diffusion on _x which pointst at mPreDensity
    SWAPPING(_x, _x0); // _x now points at mDensity
    advectDensity(_x, _x0, mU, mV); // Advection on _x which points at mDensity
}


void MyWorld::velocityStep(double *_u, double *_v, double *_u0, double *_v0) {
    SWAPPING(_u, _u0); // _u now points at mPreU
    SWAPPING(_v, _v0); // _v now points at mPreV
    diffuseVelocity(_u, _v, _u0, _v0);
    project (_u, _v, _u0, _v0);
    SWAPPING(_u, _u0); // _u now points at mU
    SWAPPING(_v, _v0); // _u now points at mV
    advectVelocity(_u, _v, _u0, _v0);
    project(_u, _v, _u0, _v0);
}

void MyWorld::diffuseDensity(double *_x, double *_x0) {
    double a = mTimeStep * mDiffusionCoef * mNumCells * mNumCells;
    linearSolve(_x, _x0, a, 1 + 4 * a, mChannels);
    setBoundary(_x);
}

void MyWorld::diffuseVelocity(double *_u, double *_v, double *_u0, double *_v0) {
    double a = mTimeStep * mViscosityCoef * mNumCells * mNumCells;
    linearSolve(_u, _u0, a, 1 + 4 * a);
    linearSolve(_v, _v0, a, 1 + 4 * a);
    setVelocityBoundary(_u, _v);
}

void MyWorld::advect(double *_d, double *_d0, double *_u, double *_v, int channels) {
    double dt0 = mTimeStep * mNumCells;  // h / x 
    for (int i = 1; i <= mNumCells; i++) {
        for (int j = 1; j <= mNumCells; j++) {
          //          if (abs(_u[IX(i, j)]) > 0.001)
            
            double x = i- dt0 * _u[IX(i,j)];  // dt0 * _u[IX(i,j)] computes how many cells can a particle travel in one time step 
            double y = j - dt0 * _v[IX(i,j)];
            if (x < 0.5) 
                x = 0.5f; 
            if (x > mNumCells + 0.5) 
                x = mNumCells + 0.5;
            int i0 = (int)x;
            int i1 = i0 + 1;
            if (y < 0.5) 
                y = 0.5;
            if (y > mNumCells + 0.5)
                y = mNumCells + 0.5;
            int j0 = (int)y;
            int j1 = j0 + 1;
            double s1 = x - i0;
            double s0 = 1 - s1;
            double t1 = y - j0;
            double t0 = 1 - t1;
            if(channels == 1) {
                _d[IX(i,j)] = s0 * (t0 * _d0[IX(i0, j0)] + t1 * _d0[IX(i0,j1)])+ s1 * (t0 * _d0[IX(i1, j0)] + t1 * _d0[IX(i1,j1)]);
            } else {
                for (int c = 0; c < channels; c++) {
                    _d[IXColor(i, j, c)] = s0 * (t0 * _d0[IXColor(i0, j0, c)] + t1 * _d0[IXColor(i0, j1, c)]) + 
                                           s1 * (t0 * _d0[IXColor(i1, j0, c)] + t1 * _d0[IXColor(i1, j1, c)]);
            }
            }
	    }
    }
}

void MyWorld::advectDensity(double *_d, double *_d0, double *_u, double *_v) {
    advect(_d, _d0, _u, _v, mChannels);
    setBoundary(_d);
}

void MyWorld::advectVelocity(double *_u, double *_v, double *_u0, double *_v0) {
    // TODO: Add velocity advection code here
    advect(_u, _u0, _u0, _v0);
    advect(_v, _v0, _u0, _v0);
    setVelocityBoundary(_u, _v);
}

void MyWorld::project(double *_u, double *_v, double *_u0, double *_v0) {
   // TODO: Add projection code here
   double h = 1.0 / mNumCells;

    for (int i = 1; i <= mNumCells; i++) {
        for (int j = 1; j <= mNumCells; j++) {
            // div
            _v0[IX(i, j)] = -0.5 * h * (_u[IX(i + 1, j)] - _u[IX(i - 1, j)] +
                                        _v[IX(i, j + 1)] - _v[IX(i, j - 1)]);
            // p
            _u0[IX(i, j)] = 0.0;
        }
    }

    //setBoundary(_v0);
    //setBoundary(_u0);
    setVelocityBoundary(_u0, _v0);

    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= mNumCells; i++) {
            for (int j = 1; j <= mNumCells; j++) {
                _u0[IX(i, j)] = (_v0[IX(i, j)] + _u0[IX(i - 1, j)] + _u0[IX(i + 1, j)] + 
                                                 _u0[IX(i, j - 1)] + _u0[IX(i, j + 1)]) / 4;
            }
        }
        setBoundary(_u0);
    }
    //setBoundary(_v);
    //setBoundary(_u);
    setVelocityBoundary(_u, _v);

    for (int i = 1; i <= mNumCells; i++) {
        for (int j = 1; j <= mNumCells; j++) {
            _u[IX(i, j)] -= 0.5 * (_u0[IX(i + 1, j)] - _u0[IX(i - 1, j)]) / h;
            _v[IX(i, j)] -= 0.5 * (_u0[IX(i, j + 1)] - _u0[IX(i, j - 1)]) / h;
        }
    }

}

void MyWorld::externalForces() {
    int size = (mNumCells + 2) * (mNumCells + 2);
    for (int i = 0; i< size; i++) {
        mPreU[i] = 0;
        mPreV[i] = 0;
    }

    std::fill_n(mPreDensity, size * mChannels, 0);
}

void MyWorld::linearSolve(double *_x, double *_x0, double _a, double _c, int channels) {
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= mNumCells; i++) {
            for (int j = 1; j <= mNumCells; j++) {
                if(channels == 1) {
                    _x[IX(i, j)] = (_x0[IX(i, j)] + _a * (_x[IX(i-1, j)] + _x[IX(i+1, j)] + _x[IX(i, j-1)] + _x[IX(i, j+1)])) / _c;
                } else {
                    for (int c = 0; c < channels; c++) {
                        _x[IXColor(i, j, c)] = (_x0[IXColor(i, j, c)] + _a * (_x[IXColor(i-1, j, c)] + 
                                                _x[IXColor(i+1, j, c)] + _x[IXColor(i, j-1, c)] + _x[IXColor(i, j+1, c)])) / _c;
                    }
                }
            }
        }
    }
}

void MyWorld::setBoundary(double *_x) {
    for (int c = 0; c < mChannels; c++) {
        for (int i = 1; i <= mNumCells; i++) {
            _x[IXColor(0 ,i, c)] = _x[IXColor(1,i, c)];
            _x[IXColor(mNumCells+1, i, c)] = _x[IXColor(mNumCells, i, c)];
            _x[IXColor(i, 0, c)] = _x[IXColor(i, 1, c)];
            _x[IXColor(i, mNumCells+1, c)] = _x[IXColor(i, mNumCells, c)];
    
        }
        _x[IXColor(0, 0, c)] = 0.5 * (_x[IXColor(1, 0, c)] + _x[IXColor(0, 1, c)]);
        _x[IXColor(0, mNumCells+1, c)] = 0.5 * (_x[IXColor(1, mNumCells+1, c)] + _x[IXColor(0, mNumCells, c)]);
        _x[IXColor(mNumCells+1, 0, c)] = 0.5 * (_x[IXColor(mNumCells, 0, c)] + _x[IXColor(mNumCells+1, 1, c)]);
        _x[IXColor(mNumCells+1, mNumCells+1, c)] = 0.5 * (_x[IXColor(mNumCells, mNumCells+1, c)] + _x[IXColor(mNumCells+1, mNumCells, c)]);
    }
}

void MyWorld::setVelocityBoundary(double *_u, double *_v) {
    for (int i = 1; i <= mNumCells; i++) {
        _u[IX(0 ,i)] = -_u[IX(1,i)];
        _u[IX(mNumCells+1, i)] = -_u[IX(mNumCells, i)];
        _u[IX(i, 0)] = _u[IX(i, 1)];
        _u[IX(i, mNumCells+1)] = _u[IX(i, mNumCells)];

        _v[IX(0 ,i)] = _v[IX(1,i)];
        _v[IX(mNumCells+1, i)] = _v[IX(mNumCells, i)];
        _v[IX(i, 0)] = -_v[IX(i, 1)];
        _v[IX(i, mNumCells+1)] = -_v[IX(i, mNumCells)];
    }
    _u[IX(0, 0)] = 0.5 * (_u[IX(1, 0)] + _u[IX(0, 1)]);
    _u[IX(0, mNumCells+1)] = 0.5 * (_u[IX(1, mNumCells+1)] + _u[IX(0, mNumCells)]);
    _u[IX(mNumCells+1, 0)] = 0.5 * (_u[IX(mNumCells, 0)] + _u[IX(mNumCells+1, 1)]);
    _u[IX(mNumCells+1, mNumCells+1)] = 0.5 * (_u[IX(mNumCells, mNumCells+1)] + _u[IX(mNumCells+1, mNumCells)]);
    _v[IX(0, 0)] = 0.5 * (_v[IX(1, 0)] + _v[IX(0, 1)]);
    _v[IX(0, mNumCells+1)] = 0.5 * (_v[IX(1, mNumCells+1)] + _v[IX(0, mNumCells)]);
    _v[IX(mNumCells+1, 0)] = 0.5 * (_v[IX(mNumCells, 0)] + _v[IX(mNumCells+1, 1)]);
    _v[IX(mNumCells+1, mNumCells+1)] = 0.5 * (_v[IX(mNumCells, mNumCells+1)] + _v[IX(mNumCells+1, mNumCells)]);

}
