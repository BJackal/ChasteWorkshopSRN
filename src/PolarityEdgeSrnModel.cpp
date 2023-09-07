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

#include "PolarityEdgeSrnModel.hpp"

PolarityEdgeSrnModel::PolarityEdgeSrnModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver)
    : AbstractOdeSrnModel(8, pOdeSolver)
{
    if (mpOdeSolver == boost::shared_ptr<AbstractCellCycleModelOdeSolver>())
    {
#ifdef CHASTE_CVODE
        mpOdeSolver = CellCycleModelOdeSolver<PolarityEdgeSrnModel, CvodeAdaptor>::Instance();
        mpOdeSolver->Initialise();
        mpOdeSolver->SetMaxSteps(10000);
#else
        mpOdeSolver = CellCycleModelOdeSolver<PolarityEdgeSrnModel, RungeKutta4IvpOdeSolver>::Instance();
        mpOdeSolver->Initialise();
        SetDt(0.001);
#endif //CHASTE_CVODE
    }
    assert(mpOdeSolver->IsSetUp());
}

PolarityEdgeSrnModel::PolarityEdgeSrnModel(const PolarityEdgeSrnModel& rModel)
    : AbstractOdeSrnModel(rModel)
{
    /*
     * Set each member variable of the new SRN model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new SRN model's member variables
     * will already have been correctly initialized in its constructor.
     *
     * Note 2: one or more of the new SRN model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new SRN model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     */
    assert(rModel.GetOdeSystem());
    AbstractOdeSystem* p_parent_system(rModel.GetOdeSystem());
    SetOdeSystem(new PolarityEdgeOdeSystem(p_parent_system->rGetStateVariables()));
    for (unsigned int i=0; i < p_parent_system->GetNumberOfParameters(); ++i)
        mpOdeSystem->SetParameter(i, p_parent_system->GetParameter(i));
}

AbstractSrnModel* PolarityEdgeSrnModel::CreateSrnModel()
{
    return new PolarityEdgeSrnModel(*this);
}

void PolarityEdgeSrnModel::SimulateToCurrentTime()
{
    // Update information before running simulation
    UpdatePolarity();
    // Run the ODE simulation as needed
    AbstractOdeSrnModel::SimulateToCurrentTime();
}

void PolarityEdgeSrnModel::Initialise()
{
    AbstractOdeSrnModel::Initialise(new PolarityEdgeOdeSystem);
}

void PolarityEdgeSrnModel::InitialiseDaughterCell()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    //A new edge is initialised with zero concentrations
    mpOdeSystem->SetStateVariable("BoundA",0.0);
    mpOdeSystem->SetStateVariable("A",0.0);
    mpOdeSystem->SetStateVariable("B",0.0);
    mpOdeSystem->SetStateVariable("C",0.0);
    mpOdeSystem->SetStateVariable("BA",0.0);
    mpOdeSystem->SetStateVariable("AB",0.0);
    mpOdeSystem->SetStateVariable("CA",0.0);
    mpOdeSystem->SetStateVariable("AC",0.0);


    mpOdeSystem->SetParameter("neighbour A",0.0);
    mpOdeSystem->SetParameter("neighbour boundA",0.0);
    mpOdeSystem->SetParameter("neighbour B",0.0);
    mpOdeSystem->SetParameter("neighbour C",0.0);
    mpOdeSystem->SetParameter("neighbour BA",0.0);
    mpOdeSystem->SetParameter("neighbour AB",0.0);
    mpOdeSystem->SetParameter("neighbour CA",0.0);
    mpOdeSystem->SetParameter("neighbour AC",0.0);

}

void PolarityEdgeSrnModel::UpdatePolarity()
{
    assert(mpOdeSystem != nullptr);
    assert(mpCell != nullptr);

    double neigh_A
    = mpCell->GetCellEdgeData()->GetItem("neighbour A")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour A", neigh_A);
    
    double neigh_boundA
    = mpCell->GetCellEdgeData()->GetItem("neighbour boundA")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour boundA", neigh_boundA);

    double neigh_B
    = mpCell->GetCellEdgeData()->GetItem("neighbour B")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour B", neigh_B);

    double neigh_C
    = mpCell->GetCellEdgeData()->GetItem("neighbour C")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour C", neigh_C);

    double neigh_BA
    = mpCell->GetCellEdgeData()->GetItem("neighbour BA")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour BA", neigh_BA);

    double neigh_AB
    = mpCell->GetCellEdgeData()->GetItem("neighbour AB")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour AB", neigh_AB);

    double neigh_CA
    = mpCell->GetCellEdgeData()->GetItem("neighbour CA")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour CA", neigh_CA);

    double neigh_AC
    = mpCell->GetCellEdgeData()->GetItem("neighbour AC")[this->GetEdgeLocalIndex()];
    mpOdeSystem->SetParameter("neighbour AC", neigh_AC);

}

double PolarityEdgeSrnModel::GetA()
{
    assert(mpOdeSystem != nullptr);
    double boundA = mpOdeSystem->rGetStateVariables()[0];
    return boundA;
}

void PolarityEdgeSrnModel::SetA(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[0] = value;
}

double PolarityEdgeSrnModel::GetBoundA()
{
    assert(mpOdeSystem != nullptr);
    double A = mpOdeSystem->rGetStateVariables()[1];
    return A;
}

void PolarityEdgeSrnModel::SetBoundA(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[1] = value;
}

double PolarityEdgeSrnModel::GetB()
{
    assert(mpOdeSystem != nullptr);
    double B = mpOdeSystem->rGetStateVariables()[2];
    return B;
}

void PolarityEdgeSrnModel::SetB(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[2] = value;
}

double PolarityEdgeSrnModel::GetC()
{
    assert(mpOdeSystem != nullptr);
    double C = mpOdeSystem->rGetStateVariables()[3];
    return C;
}

void PolarityEdgeSrnModel::SetC(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[3] = value;
}

double PolarityEdgeSrnModel::GetBA()
{
    assert(mpOdeSystem != nullptr);
    double BA = mpOdeSystem->rGetStateVariables()[4];
    return BA;
}

void PolarityEdgeSrnModel::SetBA(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[4] = value;
}

double PolarityEdgeSrnModel::GetAB()
{
    assert(mpOdeSystem != nullptr);
    double AB = mpOdeSystem->rGetStateVariables()[5];
    return AB;
}

void PolarityEdgeSrnModel::SetAB(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[5] = value;
}

double PolarityEdgeSrnModel::GetCA()
{
    assert(mpOdeSystem != nullptr);
    double CA = mpOdeSystem->rGetStateVariables()[6];
    return CA;
}

void PolarityEdgeSrnModel::SetCA(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[6] = value;
}

double PolarityEdgeSrnModel::GetAC()
{
    assert(mpOdeSystem != nullptr);
    double AC = mpOdeSystem->rGetStateVariables()[7];
    return AC;
}

void PolarityEdgeSrnModel::SetAC(double value)
{
    assert(mpOdeSystem != nullptr);
    mpOdeSystem->rGetStateVariables()[7] = value;
}

double PolarityEdgeSrnModel::GetNeighbouringBoundA() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour BoundA");
}

double PolarityEdgeSrnModel::GetNeighbouringA() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour A");
}

double PolarityEdgeSrnModel::GetNeighbouringB() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour B");
}

double PolarityEdgeSrnModel::GetNeighbouringC() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour C");
}

double PolarityEdgeSrnModel::GetNeighbouringBA() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour BA");
}

double PolarityEdgeSrnModel::GetNeighbouringAB() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour AB");
}

double PolarityEdgeSrnModel::GetNeighbouringCA() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour CA");
}

double PolarityEdgeSrnModel::GetNeighbouringAC() const
{
    assert(mpOdeSystem != nullptr);
    return mpOdeSystem->GetParameter("neighbour AC");
}


void PolarityEdgeSrnModel::OutputSrnModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractOdeSrnModel::OutputSrnModelParameters(rParamsFile);
}

void PolarityEdgeSrnModel::AddSrnQuantities(AbstractSrnModel *p_other_srn,
                                              const double scale)
{
    auto other_srn
    = static_cast<PolarityEdgeSrnModel*>(p_other_srn);
    const double other_BoundA = other_srn->GetBoundA();
    const double other_A = other_srn->GetA();
    const double other_B = other_srn->GetB();
    const double other_C = other_srn->GetC();
    const double other_BA = other_srn->GetBA();
    const double other_AB = other_srn->GetAB();
    const double other_CA = other_srn->GetCA();
    const double other_AC = other_srn->GetAC();

    const double this_BoundA = GetBoundA();
    const double this_A = GetA();
    const double this_B = GetB();
    const double this_C = GetC();
    const double this_BA = GetBA();
    const double this_AB = GetAB();
    const double this_CA = GetCA();
    const double this_AC = GetAC();

    SetBoundA(this_BoundA+scale*other_BoundA);
    SetA(this_A+scale*other_A);
    SetB(this_B+scale*other_B);
    SetC(this_C+scale*other_C);
    SetBA(this_BA+scale*other_BA);
    SetAB(this_AB+scale*other_AB);
    SetCA(this_CA+scale*other_CA);
    SetAC(this_AC+scale*other_AC);
}

void PolarityEdgeSrnModel::AddShrunkEdgeSrn(AbstractSrnModel *p_shrunk_edge_srn)
{
    // Here we assume that one half of srn quantities are endocytosed and the remaining
    // half are split between neighbouring junctions. Hence we add 1/4 of srn variables
    AddSrnQuantities(p_shrunk_edge_srn, 1.00);
}

void PolarityEdgeSrnModel::AddMergedEdgeSrn(AbstractSrnModel* p_merged_edge_srn)
{
    // Add all srn variables to this edge srn
    AddSrnQuantities(p_merged_edge_srn);
}

void PolarityEdgeSrnModel::SplitEdgeSrn(const double relative_position)
{
    //Edges with longer relative lengths after split have higher concentration
    ScaleSrnVariables(relative_position);
}


// Declare identifier for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PolarityEdgeSrnModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(PolarityEdgeSrnModel)
