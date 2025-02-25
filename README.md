# FOWLTY
FOWLTY is a MATLAB/Simulink toolbox designed to simulate offshore wind farms and introduce faults in different turbine subcomponents.

Note: The toolbox will be updated gradually, with changes documented in the Releases section.

For questions, issues, or feedback, users are encouraged to post in the Discussions section of this GitHub repository or email "ypena@mondragon.edu". However, please note that direct support is not guaranteed, and response times may vary depending on workload.

# INSTALLING THE TOOLBOX
To install and set up the FOWLTY toolbox, follow these steps:

1.- Download the latest version of the toolbox, and add the folder (and subfolders) to the Matlab path. Note that such folders (and subfolders) have to be included in the Matlab path to use FOWLTY, so it is recommended to do it every time Matlab is opened (see "AUTOMATICALLY ADD TO PATH" section at the end of the readme).

2.- From the folder, run the "install.m" file. After this, restart Matlab. There should be a "FOWLTY" section on the Simulink library browser.

# UPDATING THE TOOLBOX
Since FOWLTY is continuously being upgraded, it is recommended to download the latest files and replace the old ones. However, note that it is not necessary to run the "install.m" file again.

# USING THE TOOLBOX
To use FOWLTY, you first need to create a wind farm. To do this, open a Simulink file and go to the library browser. In the "FOWLTY" section, double-click on the "Farm Template" block and follow the instructions to generate the desired farm.

Once all the steps are completed, five different files should be generated in the selected folder:

  - a ".slx" file with the Simulink model of the farm;

  - a ".mat" file containing the wind information;
     
  - an initialisation file (".m") to initialise the required variables and run the Simulink farm model from MATLAB and visualise results (note that the model can also be run from Simulink once the initialisation file is run);
      
  - a function called "windRealisationGenerator.m" which can be used to generate more ".mat" wind files as specified in the next section;
      
  - and, if the faulty turbine model is selected, a file termed "faultScenarios.m" will also be included to define the desired faults.

Once the wind farm model is generated, it is initially linked to a single wind file corresponding to a specific mean wind speed. However, it is possible to simulate additional wind conditions by generating new wind files. This can be done using the windRealisationGenerator.m function, which allows users to specify the original wind file along with the desired parameters for the new wind realization. It is important to select the original wind file to maintain consistency, as certain wind field variables must remain unchanged for the model to function correctly.

FOWLTY also allows users to define new turbine models by modifying the nrelvals.m file located in the nrel5mw/ folder. To do so, the user must create a copy of this file and modify the relevant parameters according to the specifications of the new turbine. Once the new turbine is defined, it must be integrated into the Simulink environment by opening lib_FOWLTY.slx, duplicating an existing turbine block, and making the necessary adjustments. After running the initialization script, the user will be prompted to select the newly created turbine block in Simulink. Once confirmed, the turbine will be registered, making it available as an option when creating wind farm models.

For a comprehensive definition of the modeling assumptions and implementation details, users are referred to [1].

# REFERENCING THE WORK
When using the toolbox, please, reference [1] and [2] to acknowledge the work carried out by the authors.

[1] Peña-Sanchez, Y., Penalba, M., Knudsen, T., Nava, V., & Pardo, D. (2024). Development and validation of a health-aware floating offshore wind farm simulation platform: FOWLTY. Wind Energy and Engineering Research, vol. 2, p. 100008. https://doi.org/10.1016/j.weer.2024.100008

# ACKNOWLEDGEMENTS
Special thanks to the original developers of the "SimWindFarm" toolbox, on which the FOWLTY toolbox is based:

[2] Grunnet, J. D., Soltani, M., Knudsen, T., Kragelund, M. N., & Bak, T. (2010). Aeolus toolbox for dynamics wind farm model, simulation and control. In European wind energy conference and exhibition, EWEC 2010: Conference proceedings.

This project has received funding from the European Union’s Horizon 2020 research and innovation program under the Maire Sklodowska-Curie grant agreement No. 101034297.

# AUTOMATICALLY ADD TO PATH
One way to force MATLAB to add FOWLTY files to the path every time it opens is to create a "startup.m" file in the "MATLAB" folder located in "Documents". In that "startup.m" file, write the following line:

addpath(genpath('D:\.........\FOWLTY'))

Make sure to replace 'D:\.........\FOWLTY' with the actual location of the FOWLTY folder on your computer.
