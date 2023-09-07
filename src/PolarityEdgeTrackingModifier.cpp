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

#include "PolarityEdgeTrackingModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellSrnModel.hpp"
#include "PolarityEdgeSrnModel.hpp"

template<unsigned DIM>
PolarityEdgeTrackingModifier<DIM>::PolarityEdgeTrackingModifier()
        : AbstractCellBasedSimulationModifier<DIM>(),
        mUnboundProteinDiffusionCoefficient(0.03)
{
}

template<unsigned DIM>
PolarityEdgeTrackingModifier<DIM>::~PolarityEdgeTrackingModifier()
{
}

template<unsigned DIM>
void PolarityEdgeTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Update the cell
    this->UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarityEdgeTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarityEdgeTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Recovers each cell's edge levels proteins, and those of its neighbor's
    // Then saves them

    double D = mUnboundProteinDiffusionCoefficient;

    // For ease, store a static cast of the vertex-based cell population
    assert(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation));
    auto p_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);



    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        assert(dynamic_cast<CellSrnModel*>(cell_iter->GetSrnModel()));
        auto p_cell_srn = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());
        unsigned num_edges = p_cell_srn->GetNumEdgeSrn();

        /* Cells edge data */
        std::vector<double> BoundA_old;
        std::vector<double> A_old;
        std::vector<double> B_old;
        std::vector<double> C_old;
        std::vector<double> BA_old;
        std::vector<double> AB_old;
        std::vector<double> CA_old;
        std::vector<double> AC_old;
        std::vector<double> edge_lengths;

         for (unsigned edge_index = 0 ; edge_index  < num_edges; ++edge_index)
        {  
            // Store the current unbound protein concentrations on this edge
            auto p_edge_srn = boost::static_pointer_cast<PolarityEdgeSrnModel>(p_cell_srn->GetEdgeSrn(edge_index));
            BoundA_old.push_back(p_edge_srn->GetBoundA());
            A_old.push_back(p_edge_srn->GetA());
            B_old.push_back(p_edge_srn->GetB());
            C_old.push_back(p_edge_srn->GetC());
            BA_old.push_back(p_edge_srn->GetBA());
            AB_old.push_back(p_edge_srn->GetAB());
            CA_old.push_back(p_edge_srn->GetCA());
            AC_old.push_back(p_edge_srn->GetAC());

            // Store this edge's length
            auto p_element = p_population->GetElementCorrespondingToCell(*cell_iter);
            double edge_length = p_element->GetEdge(edge_index)->rGetLength();
            edge_lengths.push_back(edge_length);
        }

        /*
         * Update unbound protein concentrations based on a linear diffusive 
         * flux between neighbouring edges.
         */
        std::vector<double> BoundA_new(num_edges);
        std::vector<double> A_new(num_edges);
        std::vector<double> B_new(num_edges);
        std::vector<double> C_new(num_edges);
        std::vector<double> BA_new(num_edges);
        std::vector<double> AB_new(num_edges);
        std::vector<double> CA_new(num_edges);
        std::vector<double> AC_new(num_edges);

        for (unsigned edge_index = 0 ; edge_index  < num_edges; ++edge_index)
        {
           auto p_edge_srn = boost::static_pointer_cast<PolarityEdgeSrnModel>(p_cell_srn->GetEdgeSrn(edge_index));
           unsigned prev_index = (edge_index == 0) ? num_edges - 1 : edge_index - 1;
           unsigned next_index = (edge_index == num_edges - 1) ? 0 : edge_index + 1;

           ///\todo consider validity of diffusive flux expression
           //double dx = 0.5 * (edge_lengths[prev_index] + edge_lengths[edge_index]);
           double dt = SimulationTime::Instance()->GetTimeStep();
           double time_elapsed = SimulationTime::Instance()->GetTimeStepsElapsed();

            /* This if statement must be used as UpdateCellData is called for setup solve
             and at the end of each time step. Thus if no time has  elapsed diffusion should
             not be occuring */
           if (time_elapsed == 0){
           A_new[edge_index] = A_old[edge_index];
           BoundA_new[edge_index] = BoundA_old[edge_index];
           B_new[edge_index] = B_old[edge_index];
           C_new[edge_index] = C_old[edge_index];
           BA_new[edge_index] = BA_old[edge_index];
           AB_new[edge_index] = AB_old[edge_index];
           CA_new[edge_index] = CA_old[edge_index];
           AC_new[edge_index] = AC_old[edge_index];
           } else {
           A_new[edge_index] = A_old[edge_index] + D*(A_old[prev_index] - 2.0*A_old[edge_index] + A_old[next_index])*dt;
           BoundA_new[edge_index] = BoundA_old[edge_index];
           B_new[edge_index] = B_old[edge_index] + D*(B_old[prev_index] - 2.0*B_old[edge_index] + B_old[next_index])*dt;
           C_new[edge_index] = C_old[edge_index] + D*(C_old[prev_index] - 2.0*C_old[edge_index] + C_old[next_index])*dt;
           BA_new[edge_index] = BA_old[edge_index];
           AB_new[edge_index] = AB_old[edge_index];
           CA_new[edge_index] = CA_old[edge_index];
           AC_new[edge_index] = AC_old[edge_index];
           }

           p_edge_srn->SetBoundA(BoundA_new[edge_index]);
           p_edge_srn->SetA(A_new[edge_index]);
           p_edge_srn->SetB(B_new[edge_index]);
           p_edge_srn->SetC(C_new[edge_index]);
           p_edge_srn->SetBA(BA_new[edge_index]);
           p_edge_srn->SetAB(AB_new[edge_index]);
           p_edge_srn->SetCA(CA_new[edge_index]);
           p_edge_srn->SetAC(AC_new[edge_index]);

        }

        // Note: state variables must be in the same order as in PolarityOdeSystem
        cell_iter->GetCellEdgeData()->SetItem("edge boundA", BoundA_new);
        cell_iter->GetCellEdgeData()->SetItem("edge A", A_new);
        cell_iter->GetCellEdgeData()->SetItem("edge B", B_new);
        cell_iter->GetCellEdgeData()->SetItem("edge C", C_new);
        cell_iter->GetCellEdgeData()->SetItem("edge BA", BA_new);
        cell_iter->GetCellEdgeData()->SetItem("edge AB", AB_new);
        cell_iter->GetCellEdgeData()->SetItem("edge CA", CA_new);
        cell_iter->GetCellEdgeData()->SetItem("edge AC", AC_new);

    }

    //After the edge data is filled, fill the edge neighbour data
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
    {
        auto p_cell_srn = static_cast<CellSrnModel*>(cell_iter->GetSrnModel());
        unsigned num_edges = p_cell_srn->GetNumEdgeSrn();

        std::vector<double> neigh_mean_A(num_edges);
        std::vector<double> neigh_mean_BoundA(num_edges);
        std::vector<double> neigh_mean_B(num_edges);
        std::vector<double> neigh_mean_C(num_edges);
        std::vector<double> neigh_mean_BA(num_edges);
        std::vector<double> neigh_mean_AB(num_edges);
        std::vector<double> neigh_mean_CA(num_edges);
        std::vector<double> neigh_mean_AC(num_edges);

        
        for (unsigned edge_index = 0; edge_index < num_edges; ++edge_index)
        {
          auto elem_neighbours = p_population->GetNeighbouringEdgeIndices(*cell_iter, edge_index);
          for (auto neighbour : elem_neighbours)
            {
                auto p_cell = p_population->GetCellUsingLocationIndex(neighbour.first);
                auto p_data = p_cell->GetCellEdgeData();

                std::vector<double> neighbour_A_vec = p_data->GetItem("edge A");
                std::vector<double> neighbour_BoundA_vec = p_data->GetItem("edge boundA");
                std::vector<double> neighbour_B_vec = p_data->GetItem("edge B");
                std::vector<double> neighbour_C_vec = p_data->GetItem("edge C");
                std::vector<double> neighbour_BA_vec = p_data->GetItem("edge BA");
                std::vector<double> neighbour_AB_vec = p_data->GetItem("edge AB");
                std::vector<double> neighbour_CA_vec = p_data->GetItem("edge CA");
                std::vector<double> neighbour_AC_vec = p_data->GetItem("edge AC");
                
                neigh_mean_A[edge_index] += neighbour_A_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_BoundA[edge_index] += neighbour_BoundA_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_B[edge_index] += neighbour_B_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_C[edge_index] += neighbour_C_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_BA[edge_index] += neighbour_BA_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_AB[edge_index] += neighbour_AB_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_CA[edge_index] += neighbour_CA_vec[neighbour.second] / elem_neighbours.size();
                neigh_mean_AC[edge_index] += neighbour_AC_vec[neighbour.second] / elem_neighbours.size();
            }  
        }
         
        cell_iter->GetCellEdgeData()->SetItem("neighbour A", neigh_mean_A);
        cell_iter->GetCellEdgeData()->SetItem("neighbour boundA", neigh_mean_BoundA);
        cell_iter->GetCellEdgeData()->SetItem("neighbour B", neigh_mean_B);
        cell_iter->GetCellEdgeData()->SetItem("neighbour C", neigh_mean_C);
        cell_iter->GetCellEdgeData()->SetItem("neighbour BA", neigh_mean_BA);
        cell_iter->GetCellEdgeData()->SetItem("neighbour AB", neigh_mean_AB);
        cell_iter->GetCellEdgeData()->SetItem("neighbour CA", neigh_mean_CA);
        cell_iter->GetCellEdgeData()->SetItem("neighbour AC", neigh_mean_AC);
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
        
        cell_iter->GetCellEdgeData()->SetItem("in boundA", cell_iter->GetCellEdgeData()->GetItem("edge boundA"));
        cell_iter->GetCellEdgeData()->SetItem("in A", cell_iter->GetCellEdgeData()->GetItem("edge A"));
        cell_iter->GetCellEdgeData()->SetItem("in B", cell_iter->GetCellEdgeData()->GetItem("edge B"));
        cell_iter->GetCellEdgeData()->SetItem("in C", cell_iter->GetCellEdgeData()->GetItem("edge C"));
        cell_iter->GetCellEdgeData()->SetItem("in BA", cell_iter->GetCellEdgeData()->GetItem("edge BA"));
        cell_iter->GetCellEdgeData()->SetItem("in AB", cell_iter->GetCellEdgeData()->GetItem("edge AB"));
        cell_iter->GetCellEdgeData()->SetItem("in CA", cell_iter->GetCellEdgeData()->GetItem("edge CA"));
        cell_iter->GetCellEdgeData()->SetItem("in AC", cell_iter->GetCellEdgeData()->GetItem("edge AC"));

    }

}

template<unsigned DIM>
void PolarityEdgeTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PolarityEdgeTrackingModifier<1>;
template class PolarityEdgeTrackingModifier<2>;
template class PolarityEdgeTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarityEdgeTrackingModifier)
