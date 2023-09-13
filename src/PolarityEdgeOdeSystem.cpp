/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "CellwiseOdeSystemInformation.hpp"
#include "PolarityEdgeOdeSystem.hpp"

PolarityEdgeOdeSystem::PolarityEdgeOdeSystem(std::vector<double> stateVariables)
    : AbstractOdeSystem(8)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<PolarityEdgeOdeSystem>);

    /**
     * The state variables are as follows:
     *
     * 0 - Bound A homodimer concentration for this cell edge
     * 1 - Unbound A concentration for this cell edge
     * 2 - B concentration for this cell edge
     * 2 - C concentration for this cell edge
     * 4 - B-A concentration for this cell edge
     * 5 - A-B concentration for this cell edge
     * 6 - C-A concentration for this cell edge
     * 7 - A-C concentration for this cell edge
     * 
     * We store the last state variable so that it can be written
     * to file at each time step alongside the others, and visualized.
     */

    SetDefaultInitialCondition(0, 1.0); // soon overwritten
    SetDefaultInitialCondition(1, 1.0); // soon overwritten
    SetDefaultInitialCondition(2, 1.0); // soon overwritten
    SetDefaultInitialCondition(3, 1.0); // soon overwritten
    SetDefaultInitialCondition(4, 1.0); // soon overwritten
    SetDefaultInitialCondition(5, 1.0); // soon overwritten
    SetDefaultInitialCondition(6, 1.0); // soon overwritten
    SetDefaultInitialCondition(7, 1.0); // soon overwritten

    this->mParameters.push_back(0.5);
    //By default zero.
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);
    this->mParameters.push_back(0.0);
    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

PolarityEdgeOdeSystem::~PolarityEdgeOdeSystem()
{
}

void PolarityEdgeOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    const double BoundA = rY[1];
    const double A = rY[0];
    const double B = rY[2];
    const double C = rY[3];
    const double BA = rY[4];
    const double AB = rY[5];
    const double CA = rY[6];
    const double AC = rY[7];

    const double neigh_A = this->mParameters[0]; // Shorthand for "this->mParameter("neighbor delta");"
    const double neigh_B = this->mParameters[2]; // Shorthand for "this->mParameter("neighbor notch");"
    const double neigh_C = this->mParameters[3]; // Shorthand for "this->mParameter("neighbor notch");"
    const double neigh_BA = this->mParameters[4]; // Shorthand for "this->mParameter("neighbor notch");"
    const double neigh_CA = this->mParameters[6]; // Shorthand for "this->mParameter("neighbor notch");"

    const double KD1 = 5;
    const double KD2 = 0.1;

    const double k = 1.0;
    const double v1 = KD1*k;
    const double v2 = KD2*k;

    const double K = 0.1665;
    const double VF = 10.0;
    const double VS = 10.0;
    const double xF = BA;
    const double xS = CA;
    const double xFm = (neigh_BA);
    const double xSm = (neigh_CA);
    const double w = 2.0;

    const double hF = 1 + (((VF -1)*pow(xF,w))/(pow(K,w) + pow(xF,w)));
    const double hS = 1 + (((VS -1)*pow(xS,w))/(pow(K,w) + pow(xS,w)));
    const double hFm = 1 + (((VF -1)*pow(xFm,w))/(pow(K,w) + pow(xFm,w)));
    const double hSm = 1 + (((VS -1)*pow(xSm,w))/(pow(K,w) + pow(xSm,w)));

    const double R1 = k*(A*neigh_A) - (v1*(BoundA)); // A -> Bound A

    const double R2 = k*(B*BoundA) - v2*hS*xS*(BA); // FZ:FL -> B + Bound A
    const double Rm2 = k*(neigh_B*BoundA) - v2*hSm*xSm*(AB); // FL:FZm -> Bm + Bound A

    const double R3 = k*(C*BoundA) - v2*hF*xF*(CA); // Stb:FL -> B + Bound A
    const double Rm3 = k*((neigh_C)*BoundA) -v2*hFm*xFm*(AC);// FL:Stbm -> Bm + Bound A


    rDY[1] =  R1 - R2 -Rm2 -R3 -Rm3;// d[Bound A]/dt

    rDY[0] = -R1 ; // d[A]/dt
    
    rDY[2] = -R2; // d[B]/dt

    rDY[3] = -R3; // d[C]/dt

    rDY[4] = R2; // d[FZ:FL:FL]/dt

    rDY[5] = Rm2; // d[FL:FL:FZ]/dt

    rDY[6] = R3; // d[STB:FL:FL]/dt

    rDY[7] = Rm3; // d[FL:FL:STB]/dt

}

template<>
void CellwiseOdeSystemInformation<PolarityEdgeOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("BoundA");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0); // will be filled in later

    this->mVariableNames.push_back("A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0); // will be filled in later

    this->mVariableNames.push_back("B");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0); // will be filled in later

    this->mVariableNames.push_back("C");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0); // will be filled in later

    this->mVariableNames.push_back("BA");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0); // will be filled in later

    this->mVariableNames.push_back("AB");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0); // will be filled in later

    this->mVariableNames.push_back("CA");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0); // will be filled in later

    this->mVariableNames.push_back("AC");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0); // will be filled in later

    this->mParameterNames.push_back("neighbour A");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour boundA");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour B");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour C");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour BA");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour AB");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour CA");
    this->mParameterUnits.push_back("non-dim");
    this->mParameterNames.push_back("neighbour AC");
    this->mParameterUnits.push_back("non-dim");

    this->mInitialised = true;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PolarityEdgeOdeSystem)
