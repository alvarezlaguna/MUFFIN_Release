//
//  LinearSolver.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 24/09/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//

#ifndef LinearSolver_hpp
#define LinearSolver_hpp

#include <stdio.h>
#include <vector>
#include "../Component.hpp"
#include <string>


using namespace std;

typedef vector<vector<double> > Matrix;

class LinearSolver : public Component
{
public:
    LinearSolver(string name, unsigned int size) : Component(name)  {}
    virtual ~LinearSolver() {}
    void setup(){}
    virtual void solveLinearSystem(const Matrix& A, vector<double>& x, const vector<double>& B) = 0;
};


#endif /* LinearSolver_hpp */
