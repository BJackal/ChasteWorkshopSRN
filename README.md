# Modelling cell-cell signalling with state reaction networks in Chaste

So-called 'sub-cellular reaction networks' (SRNs) can be used to mathematically model subcellular biological processes such as gene regulatory networks, protein-protein interactions, and signalling pathways. 
In Chaste, these processes may influence, and/or be influenced by, cellular behaviours such as proliferation, differentiation, adhesion, and motility, as well as cell-cell communication.
In this practical, you will investigate how to simulate two examples of cell-cell communication using SRNs in Chaste, building on the previous vertex model examples you looked at on Monday and Tuesday.

First you need to follow the steps outlined in the workshop on Tuesday Chaste to correctly copy and initalise a user project: https://github.com/Chaste/IntroToChasteCpp

 * If you are running Chaste through docker you will first need to start up your image
 * Navigate to your Chaste projects folder: ```cd /path/to/Chaste/projects``` for native usage or within docker simply navigate to ```Chaste/projects``` for your docker image after building Chaste navigate to chaste/projects.
 * Clone this user project: ``` git clone https://github.com/BJackal/ChasteWorkshopSRN.git```
 * Navigate to the cloned user project from the terminal and then in the temrinal execute ```python setup_project.py``` to setup the project. This user project only relies on cell-based Chaste, so answer `Y` to cell-based and `n` to lung and heart. 
 * Now you are ready to make sure that the user project will compile correctly. To do this, use the command line to navigate to the Chaste build folder. If you are running Chaste nativley this should be made outside of the main Chaste source folder. If you are running in Docker you will need to navigate to the build folder within your Docker image after running the build script. Within the build directory at the terminal run: ```ccmake path/to/Chaste```. While this is running, you will be prompted for an input to configure your user project. Press ```c ```, then on completion press ```e``` to exit. Next, configure with ```c``` once more, followed up by ```e``` to exit and finally ```g``` to generate. You are now ready to to make our tests.
 * To ensure that this has been done correctly, you can build the example test called `TestHello` and ensure that it passes. First, while within your build directory run ```make TestHello_ChasteWorkshopSRN``` and on completion you can run the example test using ```ctest -V -R TestHello_ChasteWorkshopSRN```. If this test succesfully makes and passes, youe know that the user project is compiled correctly. 

Here you will explore two examples of using SRN models in chaste. First, you will explore an example in Chaste of the Delta-Notch signalling pathway. Here you will utilise a modified version of Collier et als (1996, https://doi.org/10.1006/jtbi.1996.0233) work for the feedback loop between Delta and Notch. The one major difference in our system being that Delta concentrations in the neighbouring cells are used directly. Second, you will be able to investigate the formation of planar polarity in a 2D tissue based on a simplified version of the work of Fisher et al (2019, https://doi.org/10.1016/j.isci.2019.06.021).

## Example 1: Notch signalling
The test you will be following for this section can be found within ```path/to/Chaste/projects/ChasteWorkshopSRN/test/TestDeltaNotchSRN```. First, open this file in visual studio code. You can look over the commented code describing what each line of text is being used for. 

All files will start out with relavent includes for which parts of Chaste you wish to utilise and this is also where you would include any of your own custom code that you may use in your own simulations.

For the first run through you will look at an edge based SRN model for a line of 70 cells. To run this simulation you first need to navigate to your ```/build``` directory.  Then in your console input ```make -j4 TestDeltaNotchSRN``` to make your simulation. Any time you alter your test you need to remember to make your test again. On completion execute ```ctest -V -R TestDeltaNotchSRN```, this will execute the test file and run your simulation. Once your test runs successfully you will need to utilise a tool such as Paraview to explore your results. 

If you do not already have Paraview installed on your system, you can download it from the Paraview website here: https://www.paraview.org/download/

In Paraview first navigate to File -> Open. 

![File](https://user-images.githubusercontent.com/44051158/267284363-fd8975ab-6d5e-4b08-a0c4-74bce976413d.png) 

Then navigate to your test results folder, by default if running Chaste nativley this should be /tmp/YourUserName/testoutput/TestDeltaNotchEdgeODESimulation. If running Chaste through docker these outputs should be under /testoutput. Opening the results_from_time_0.vtu file. Click the green apply button for your cells to appear. 

![Apply](https://user-images.githubusercontent.com/44051158/267284370-dbd9b8c6-819a-4be9-859a-e33cbaee02e5.png)

Then on the second row of the toolbar there will be two drop down menus defaulted to "Solid colour" and "Surface". You will need to change these to "edge delta" and "Surface With Edges". 

![EdgeData](https://user-images.githubusercontent.com/44051158/267284375-bd646b3d-4f8c-4358-a60f-6fa24bc2da41.png)

Finally, before playing the simulation you need to modify your colour bar for the entire time series of the simulation. To do this you can click the 6th button on the second row of the toolbars ( This should be a green arrow pointing in two directions with a small t) then click rescale when prompted withthe disable automatic rescaling option. 

![ColourBarTimeScaling](https://user-images.githubusercontent.com/44051158/267284376-13b09ccf-77d7-4a31-a20f-772338166f39.png)

Now you can play your simulation and watch how the edge based Delta evolves over time in the simulation. Similarlly, you can also change from your "edge delta" to "edge notch" to view how that also changes in time. 

If you wish to explore and plot how a given edges cocnentration of Delta and Notch changes in time you need to first click the "Select Cells On" button (This looks like a green triangle surrounded by a dotted box) then select your given edge of interest. In this example right most edge of the second cell conencting to its neigbhour has been selected. You can tell what you have selected by a pruple border around the area of interest. 

![EdgeSelection](https://user-images.githubusercontent.com/44051158/267284377-3735d07b-3256-4d63-8df5-1e14acdcacfd.png)

Now, navigate to Filters -> Data Analysis -> Plot Selection Over time. Then click the green apply button.This will create a graph showing all parameters on that edge.

![GraphFirst](https://user-images.githubusercontent.com/44051158/267284352-040fdefb-e231-40a1-b6a1-893a3d1ec95d.png)

For now turn off all parameters except edge notch and edge delta and you should have a graph that looks like the following image.

![GraphSecond](https://user-images.githubusercontent.com/44051158/267284359-fc11e6d5-e2e9-4709-ad82-7788d14c4e16.png)

Now you can see how your concentrations of Delta and Notch evolve on this edge in time. You can try this across multiple edges between neighbours to see how each individual edge evolves in time.

For the second run of the simulation you can change your cells to be in a small tissue instead of a single line. To do this you will need to change line 101 in your test file from ```HoneycombVertexMeshGenerator generator (70,1)``` to  ```HoneycombVertexMeshGenerator generator (6,6)```. This will construct your cells in a 6 by 6 grid of cells instead of a single line. Re-run the simulation again using ```make -j4 TestDeltaNotchSRN``` and  ```ctest -V -R TestDeltaNotchSRN``` then view how this looks in paraview.  

In some case you may want to add proliferation to your cells. To do this you will need to change a few sections of the underlying test file. First, you will need to change line 114 from ```NoCellCycleModel* p_cc_model = new NoCellCycleModel``` to ```UnifromG1GenerationCellCycleModel* p_cc_model = new UnifromG1GenerationCellCycleModel```. This means your cells will now have a stochastic cell cycle model, allowing for your cells to grow and divide based on their current phase. This will assign your cells phases for growth and division events from a uniform distribution. In addition, you will also need to add a target area modifier for your cells, to do this uncomment line 195 and 196 in your TestDeltaNotchSRN file. Finally, you will also need to add a force to your simulation to give your cells some movement, to do this uncomment line 192 and 193. This will pass a Farhardifar force to your cells. This is not strictly necessary for cell dvision. However, if you do not pass some kind of force to your simulation cells will simply divide in place creating a spider web of cells within the initial cell volume.

Now you are ready to run your simulation again. Still in your ```/build``` directory execute ```make -j4 TestDeltaNotchSRN``` and  ```ctest -V -R TestDeltaNotchSRN``` again as previously. This may take a little longer to run than before due to the addition of cell proliferation in your system. Viewing the simulation results in paraview you will be able to see a much more dynamic system than in your first run case. In addition to cell movement and division you should also be able to find various types of cell transistions (T1, T2 etc ....) as some cells will exhange neighbours or be removed from the system as the tissue expands.

You can further explore this test by trying different cell configurations. Simillarly you can also investigate how changing the initial conditions in your simulation, such as Delta and Notch, may lead to different end states for your simulation.

## Example 2: Planar polarity

Another use case for edge based SRN simulations is modelling planar polarity in epethelial tissues. Planar polarity(PP) is a form of protein mediated patterning across tissue due to signalling. Creating distinctive patterns across large length scales in many tissues. For example, in the wing of the Drasophila fruit fly wings there are many hairs across the wings used for sensory pruposes. These hairs have been noted to be all pointing from the proximal side of wings out towards the distal end. It has been found that there is a set of core proteins in the Drosophila wings which which create a simillar patterning on the cellular to tissue scale. This pattern is formed by two main proteins Frizzled(Fz) and Strabismus (Stb). These proteins then localise proximally and or distally based on positive feedback with their own species and ihibit the binding of the opposite protein on the same edge. Fz and Stb binding at cellular junctions to Flamingo (Fl) homodimers which bind across junctions.

In this example, you will explore a simplified version of this model focusing on only three "protein" species: A, B and C. In this case protein A can form homodimers with its own species in neighbouring cells edges and proteins B and C can then bind to these homodimers. First, you will need to open the Polarity test. This can be found in  ```path/to/Chaste/projects/ChasteWorkshopSRN/```. In addition, as you will be utilising some code for this simulation outside of the Chaste source code you can also find several additional polarity files within ```path/to/Chaste/prohects/ChasteWorkshopSRN/src```. It is good practice to always keep any of your own code which you will be editing as part of your project within this section and not directlly editing underlying Chaste source code. This is so if any issues arise it will be easier to ascertain if the new code or something else is causing the issue.

With your Polarity test open you can see that again we begin with all of our relavent includes. In addition, all of custom code outside of the Chaste trunk that this tests utilises is also included here. The test utilises a regular vertex mesh the same as the previous example with an initial configuration of 6 by 6 cells. In this test you can see that unbound proteins are initialised at none-zero values and any bound proteins are initialised at a concentration of zero. 

Following previous steps, you can run the Polarity test by executing ```make -j4 TestPolaritySRN``` and on completion execute ```ctest -V -R TestPolaritySRN```. Then you can launch the simulation results in paraview and investigate how polarity forms across the tissue in time. You can use the previous methods described in Example 1 to investigate the evolution of the protein species over time and investiagte how the species interact with each other.

You can also edit the underlying source code for your test in the Chaste/projects/ChasteWorkShopSRN/src folder. For example, you can open PolarityEdgeTrackingModifier.cpp. This file contains several underlying functions for keeping track of cellular edge data. On line 38 you will find the class variable```mUnboundProteinDiffusionCoefficient(0.03)```. This represents the diffusion rates for your proteins and is utilised in the function ```UpdateCellData``` for diffusion of proteins around your cells edges. Try reducing or increasing this coefficient then save your changes to the file. Now that you have updated this underlying src file you will need to cmake our build folder again. To do this lets go back to our build folder (outside the main chaste soruce code folder) and run at command line ```cmake /path/to/Chaste/``` then ```make -j4 TestPolaritySRN``` and on completion execute ```ctest -V -R TestPolaritySRN```. Opening up the results in paraview what change if any do you find to the temoporal behaviour of the proteins compared to the original diffusion coefficient ?

In all of the tests so far cells have been on a uniform vertex mesh. In some cases you may wish to investigate the dynamics of cells on a more randomly distrubted mesh. For example, you could create cells which are randomly distibuted with their edges defined by a Voronoi teselation. To do this, re-open your ```TestPolaritySRN.hpp``` and comment out line 111 and repalce this with line 112. This will switch out your uniform vertex mesh for a Voronoi vertex mesh. This will create a dynamic mesh with cells edges defined based on distances from each others centres. Now re-run your simulation to view how this looks in paraview and compare it to the previous uniform mesh.

Simillar to in Example 1 you can also repalce your ```NoCellCycleModel* p_cc_model = new NoCellCycleModel``` with ```UnifromG1GenerationCellCycleModel* p_cc_model = new UnifromG1GenerationCellCycleModel```. You will also need to uncomment lines 202 - 206. This will allow your cells to grow and divide. On division, the new daughter cells edges will receive a portion of the concentration that existed on the dividing mother cells edge. By default this will be  

You could further explore this test by altering the initial conditions of your A,B and C concentrations. For example, how does your final simulation result change if you give none equal cocnentrations of the initial proteins ? How is the end result effected if you "flooded" the simulation with a higher initial cocnentraiton of the A homodimer forming protein. You could also try increasing or decreasing the initial bias caused by the offset. How does altering this effect your end results.

# Modelling cell-environment feedbacks with partial differential equations in Chaste

Reaction-advection-diffusion partial differential equations (PDEs) can be used to mathematically model paracrine signalling, nutrient transport, and other forms of cell-environment feedbacks. 
In Chaste, dynamic or steady-state PDEs (or PDE systems) can be defined on the spatial domain occupied by a population of cells. 
Simulation modifiers are used to numerically integrate these PDEs over time, and their numerical solution may be used by cells to inform their behaviour via `CellData` items. 
In this practical, you will investigate two examples of cell-environment feedbacks in Chaste, oxygen-dependent tumour spheroid growth and morphogen-dependent cell proliferation.

## Example 1: Oxygen-dependent tumour spheroid growth

Read through the [https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials/RunningTumourSpheroidSimulations](tumour spheroid 'user tutorial') that can be found on the old Chaste wiki. 
If you are able to explore your local copy of the Chaste source code, e.g. in Visual Studio code, then you can find this user tutorial as a test suite located at `cell_based/test/tutorial/TestRunningTumourSpheroidSimulationsTutorial.hpp`. 
Try running this yourself for different values of the oxygen uptake rate and the boundary oxygen concentration. 
You can also explore varying the oxygen-dependent cell-cycle model parameters by calling the methods `SetHypoxicConcentration()`, `SetQuiescentConcentration()`, or `SetCriticalHypoxicDuration()` on `p_model`.
If you're not sure how to do this after having read the user tutorial, then please ask one of the demonstrators.

## Example 2: Morphogen-dependent proliferation

Visit the [https://chaste.cs.ox.ac.uk/trac/wiki/PaperTutorials/CellBasedComparison2017/MorphogenMonolayer](morphogen-dependent proliferation 'paper tutorial') that can be found on the old Chaste wiki. 
Try running this yourself locally, by creating a new test suite and pasting the code found at the bottom of the wiki page into it. 
If you're not sure how to create and run a new test suite after having read the guidance [https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials/WritingTests](here), then please ask one of the demonstrators.
