/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTHELLO_ANOTCH_HPP_
#define TESTHELLO_ANOTCH_HPP_

#include <cxxtest/TestSuite.h>
/* Most Chaste code uses PETSc to solve linear algebra problems.  This involves starting PETSc at the beginning of a test-suite
 * and closing it at the end.  (If you never run code in parallel then it is safe to replace PetscSetupAndFinalize.hpp with FakePetscSetup.hpp)
 */
#include "AbstractCellBasedTestSuite.hpp"
#include "CellAgesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellSrnModel.hpp"
#include "CellVolumesWriter.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DeltaNotchEdgeSrnModel.hpp"
#include "DeltaNotchEdgeTrackingModifier.hpp"
#include "DeltaNotchInteriorSrnModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "FarhadifarForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"
#include "PetscSetupAndFinalize.hpp"

/**
 * @file
 *
 * This is an example of a test suite, used to test the source
 * code, and also used to run simulations for A Notch based simulations.
 * Here our A and notch will evolve on the edge of our cells based on 
 * the A-Notch ODE system described by Collier et al,
 * "Pattern formation by lateral inhibition with feedback: a mathematical
 * model of A-notch intercellular signalling" (Journal of Theoretical
 * Biology 183:429-446, 1996).
 * 
 * Here, however, we include edge based model: A and Notch interactions between each cell
 * are modelled directly. We use similar ODE system as by Collier et al., except that we modify terms
 * corresponding to means of neighbour concentrations of A/Notch.
 *
 * You can #include any of the files in the project 'src' folder.
 * For example here we #include "CellVolumeWriter.hpp"
 *
 * You can utilise any of the code in the main the Chaste trunk
 * in exactly the same way.
 * NOTE: you will have to alter the project SConscript file lines 41-44
 * to enable #including of code from the 'heart', 'cell_based' or 'crypt'
 * components of Chaste.
 */

class TestHello_DeltaNotch : public AbstractCellBasedTestSuite
{
public:
    void TestRunningMultiODECellWithEdges()
    {
        /* First we create a regular vertex mesh containg a mesh of 70 x 1 elements. */
        HoneycombVertexMeshGenerator generator(70, 1);
        boost::shared_ptr<MutableVertexMesh<2, 2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        /*The mutation and prolfierative type of our cells also needs to be declared*/
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_diff_type);
        
        /*Now we need to loop over all the elements in the mesh to set up our initial conditions*/
        for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            /* To be begin with, we will initalise our cells with the dummy no cycle 
                model to prevent proliferation*/
            NoCellCycleModel* p_cc_model = new NoCellCycleModel();
            p_cc_model->SetDimension(2);

            /*Get the current element*/
            auto p_element = p_mesh->GetElement(elem_index);

            /* Initialise edge based SRN */
            auto p_cell_edge_srn_model = new CellSrnModel();

            /* We choose to initialise the total concentrations to random levels */
            auto delta_concentration = RandomNumberGenerator::Instance()->ranf();
            auto notch_concentration = RandomNumberGenerator::Instance()->ranf();
            
            /*Our concentrations will be initalised based on the edges current length*/
            double total_edge_length = 0.0;
            for (unsigned i = 0; i < p_element->GetNumEdges(); i++)
            {
                total_edge_length += p_element->GetEdge(i)->rGetLength();
            }

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges(); i++)
            {
                /*Get the current edge of the element at this indexx*/
                auto p_elem_edge = p_element->GetEdge(i);
                /*Get the length of that edge*/
                auto p_edge_length = p_elem_edge->rGetLength();
                std::vector<double> initial_conditions;

                /* Initial concentration of delta and notch vary depending on the edge length */
                initial_conditions.push_back(p_edge_length / total_edge_length * delta_concentration);
                initial_conditions.push_back(p_edge_length / total_edge_length * notch_concentration);
                
                /*Now we can set our initial conditions and add an SRN model to the edge*/
                MAKE_PTR(DeltaNotchEdgeSrnModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }
            
            /*Now we are ready to begin setting up our cells passing their cell cycle model,
              edge SRN model and proliferative type*/
            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            
            /*We also must give our cells a birth time, here we randomly initialise it*/
            const double birth_time = -RandomNumberGenerator::Instance()->ranf() * 12.0;
            p_cell->SetBirthTime(birth_time);
            /*Finally we push back the current cell to our overall cells*/
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        /*All information we wish to be output to file must be called through our writers.
        For example we here call to add the cells volumes*/
        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. First we pass our cell population as an off lattice simulation.*/
        OffLatticeSimulation<2> simulator(cell_population);
        /*Then we tell our simulation were to output its results, by default this goes to /tmp/YourUser/testoutput*/
        simulator.SetOutputDirectory("TestDeltaNotchEdgeOnlyODESimulation1");
        /*Defining that our simulation will output its deta every tenth Dt, with a Dt of 0.001 and an end time of 100*/
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetDt(0.1);
        simulator.SetEndTime(500);

        /* We pass our edge tracking modifier to the cell to call Update CellEdgeData, this allows 
         the SRN simulation to be initalised correctlly*/
        MAKE_PTR(DeltaNotchEdgeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);
        
        /*We can also pass a force modifier to our simulation but will initially leave this turned off*/
        //MAKE_PTR(FarhadifarForce<2>, p_force);
        //simulator.AddForce(p_force);

        //MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        //simulator.AddSimulationModifier(p_growth_modifier);
        
        /*Finally, we run our simulation and assert that the solver throws nothing*/
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }
};

#endif /*TESTHELLO_ANOTCH_HPP_*/
