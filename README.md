# Delta-Notch and cell polarity SRN modelling project for use with Chaste.

Biological systems are an extensive research area across multiple disciplines. The high level of complexity within these systems due to dependencies on multiple factors makes them challenging to investigate. However, many of these systems can be represented as a network of biochemical reactions. For instance Subcellular Reaction Networks (SRNs) can be used to model multiple biological systems, such as: gene regulation networks, predicting protein-protein interaction networks and many more.

First we need to follow the steps outlined on the chaste website for correctlly initalizing our user project: 
https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/UserProjects

Thus we first must navigate to our Chaste project folder: ```cd /path/to/Chaste/projects```
Then we must clone this project: ``` git clone https://github.com/BJackal/DeltaNotch.git```
Now we navigate into our cloned user project from the command line and run our setup: ```python setup_project.py```
Here our project only relies on the cell based portion of Chaste so we answer yes to cell based and no lung and heart.
Now we are ready to make sure that our project will compile correctly. To do this we use the command line to navigate to our chaste build folder (This should be outside of the main Chaste source folder) and run: ```ccmake path/to/Chaste```
While this is running we will be prompted for an input to configure our project. For this we will press ```c ``` then on completion ```e```` to exit.
Then we will configure with ```c``` onece more, followed up by ```e``` to exit and finally ```g``` to generate.
We are no ready to to make our tests.

To ensure that this has been done correctlly we can build our example test of TestHello and ensure that it passes.
First at the command line we run ```make TestHello``` and on completion we can run the example test using ```ctest -R TestHello$```.
If this test succesfully makes and pass we know that our project is correctlly compiled and working.

Here we will explore two examples of using SRN models in chaste. First, we will explore an example in Chaste of the Delta-Notch signalling pathway. Here we utilise the ODE system presented in the work Collier et al (1996) for the feedback loop between Delta and Notch. The one major difference in our system being that Delta concentrations in the neighbouring cells are used directly.

# Example 1
The test we will be followingfor our first example can be found within /test/Test_DeltaNotch. Let's open this file in visual studio code. First we can look over the commented code describing what each line of text is being used for. All files will start out with our relavent includes for which parts of Chaste we wish to utilise and this is also where you woul include any of your own custom code that you may use in your own simulations.

For our first run through we will look at an edge based SRN model for a line of 70 cells. To run this simulation we simply need to go to our command line and input ```make -j4 Test_DeltaNotch``` and on completion execute ```ctest -V -R Test_DeltaNotch```. To explore our results we will utilise Paraview. In Paraview first navigate to File -> Open. Then navigate to our results, by default this should be /tmp/YourUserName/testoutput/TestDeltaNotchEdgeODESimulation. Opening the results_from_time_0.vtu file. Clicl the green apply button for your cells to appear. Then on the second row of the toolbar there will be two drop down menus defaulted to "Solid colour" and "Surface". Lets change these to "edge delta" and "Surface With Edges". Finally, before playing our simulation we need to modify our colour bar for the entire time series of the simulation. To do this we can click the 6th button on the second row of the toolbars ( This should be a green arrow pointing in two directions with a small t) then click rescale when prompted withthe disable automatic rescaling option. Now we can play our simulation and watch how the edge based Delta evolves over time in the simulation. Similarlly, we can also change from our "edge delta" to "edge notch" to view how that also changes in time. 

If we wished to explore and plot how a given edges cocnentration of delta and notch changed in time we need to first click the "Select Cells On" button (This looks like a green triangle surrounded by a dotted box) then select our given edge. here lets click left hand edge of the second cell conencting to its neigbhour. Then simply navigate to Filters -> Data Analysis -> Plot Selection Over time. Then click the green apply button.This will create a graph showing all parameters on that edge, for now lets turn off all parameters except edge notch and edge delta. Now we can see how our concentrations of Delta and Notch evolve on this edge in time. Try this across multiple edges between neighbours to see how each individual edge evolves in time.





If you clone this repository, you should make sure to rename the template_project folder with your project name and run the 'setup_project.py' script to avoid conflicts if you have multiple projects.
