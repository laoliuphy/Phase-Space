//
//  PhaseSpace.hpp
//  bb
//
//  Created by Yandong Liu on 3/14/16.
//  Copyright © 2016 Yandong Liu. All rights reserved.
//

#ifndef PhaseSpace_hpp
#define PhaseSpace_hpp

#include <stdio.h>
#include "Calculator.hpp"
#include "Spinor.hpp"
#include "ran.h"


void Two_Body_P(vector<Particle>&, vector<Particle>&, vector<double>&);/*2->2*/

void Two_Body_P(vector<Particle>& a, vector<Particle>& c , vector<double>& b , double m , double width);

void Two_Body_P(int, vector<Particle>&, vector<Particle>&, vector<double>&);/*a->2*/

void Three_Body_P(vector<Particle>&, vector<Particle>&, vector<double>&);/*2->3*/

void Three_Body_P(vector<Particle>&, vector<Particle>&, vector<double>&, double, double);

void Three_Body_P(int , vector<Particle>&, vector<Particle>&, vector<double>&); /*a->3*/

void Three_Body_P(int, vector<Particle>&, vector<Particle>&, vector<double>&, double mass, double width); /*a->3 with a mediator of mass m and width */


double Two_Body(vector<Particle>&, vector<Particle>&, vector<double>&); /*2->2 the weight of phase space */

double Two_Body(int, vector<Particle>&, vector<Particle>&, vector<double>&);


double Three_Body(vector<Particle> & , vector<Particle> & , vector<double> &); /*2->3 the weight of the phase space point*/

double Three_Body(int , vector<Particle> & , vector<Particle> & , vector<double> &);

double Three_Body(int, vector<Particle>&, vector<Particle>&, vector<double>&, double, double);

#endif /* PhaseSpace_hpp */
