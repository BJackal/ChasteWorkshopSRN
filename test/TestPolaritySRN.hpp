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
#include "PetscSetupAndFinalize.hpp"

#include "AbstractCellBasedTestSuite.hpp"

/* Here, we include the necessary writers here for late outputting information to programmes such as paraview. For example we use #include "CellAgesWriter.hpp"
*which allows use to get data on how long a cell has been in our simulation.*/

#include "CellAgesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellSrnModel.hpp"
#include "CellVolumesWriter.hpp"

#include "CheckpointArchiveTypes.hpp"

/*As we wish to simulate our Edge based A Notch system we need to include the relavent files here.*/
#include "PolarityEdgeSrnModel.hpp"
#include "PolarityEdgeTrackingModifier.hpp"

/*Here we include all of the relavent other files to be used throughout or simulation setup
* These will be discussed further into the file.*/
#include "HoneycombVertexMeshGenerator.hpp"
#include "VoronoiVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "FarhadifarForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "NoCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"

/**
 * @file
 *
 * This is an example of a test suite, used to test the source
 * code, and also used to run simulations for a Poalrity based simulations.
 * Here our proteins will evolve and diffuse on the edge of our cells based on 
 * a simplified version of Fisher et als work,
 * "Experimental and Theoretical Evidence for Bidirectional Singaling via Core
 Planar Polarity Protein Complexes in Drosophila" (iScience  17:49-66, 2019).
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

class TestPolaritySRN : public AbstractCellBasedTestSuite
{
public:
    void TestPolarity()
    {
        /* First we create a regular vertex mesh. */
        HoneycombVertexMeshGenerator generator(6, 6);
        // VoronoiVertexMeshGenerator generator(3, 3,1,1.0);
        boost::shared_ptr<MutableVertexMesh<2, 2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            /* Initalise cell cycle */
            NoCellCycleModel* p_cc_model = new NoCellCycleModel();
            p_cc_model->SetDimension(2);

            auto p_element = p_mesh->GetElement(elem_index);

            /* Initialise edge based SRN */
            auto p_cell_edge_srn_model = new CellSrnModel();
            /* We choose to initialise the total concentrations to random levels */
            auto A_concentration = 0.333;
            auto BoundA_concentration = 0.0;
            auto B_concentration = 0.333;
            auto C_concentration = 0.333;
            auto AB_concentration = 0;
            auto BA_concentration = 0;
            auto AC_concentration = 0;
            auto CA_concentration = 0;

            /* Gets the edges of the element and create an SRN for each edge */
            for (unsigned i = 0; i < p_element->GetNumEdges(); i ++)
            {
                std::vector<double> initial_conditions;

                double offset1 = 1.000;
                //double offset2;

                if( i==0  || i==5){
                  offset1 = 0.999;
                }else if(i==1 || i==4){
                  offset1 = 1.000;
                }else if( i== 2 || i == 3){
                  offset1 = 1.001;
                }
               
                initial_conditions.push_back(A_concentration);
                initial_conditions.push_back(BoundA_concentration);
                initial_conditions.push_back((B_concentration * offset1));
                initial_conditions.push_back(C_concentration);
                initial_conditions.push_back(BA_concentration);
                initial_conditions.push_back(AB_concentration);
                initial_conditions.push_back(CA_concentration);
                initial_conditions.push_back(AC_concentration);
                }
                MAKE_PTR(PolarityEdgeSrnModel, p_srn_model);
                p_srn_model->SetInitialConditions(initial_conditions);
                p_cell_edge_srn_model->AddEdgeSrnModel(p_srn_model);
            }

            CellPtr p_cell(new Cell(p_state, p_cc_model, p_cell_edge_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);


            double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        // cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();
        // cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPolarityEdgeOnlyODESimulation1");
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetDt(0.1);
        simulator.SetEndTime(100);

        /* Update CellEdgeData so that SRN simulations can run properly */
        MAKE_PTR(PolarityEdgeTrackingModifier<2>, p_modifier);
        
        simulator.AddSimulationModifier(p_modifier);

        //MAKE_PTR(FarhadifarForce<2>, p_force);
        //simulator.AddForce(p_force);

        //MAKE_PTR(TargetAreaLinearGrowthModifier<2>, p_growth_modifier);
        //simulator.AddSimulationModifier(p_growth_modifier);

        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }
};

#endif /*TESTHELLO_ANOTCH_HPP_*/
