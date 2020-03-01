# ShockDetecionMatlabCode

Code samples that have been used in the Master Project for importing the data to the PC. Almost all below Matlab scripts have been used to get the data from different hardware (microcontroller, Oscilloscope) and display the data in a plot for evaluating the measurement, Except the first one:

ComparisonOfLegacyData: Opens the two mat. files, where the acceleration measurement is stored, trim the part with the noise, stores the values in txt format for processing in different programming environments, calculates the energy of the impulse and stores the plot in a jpeg. This script was used to look into the measurement data from the previous system that has been developed from the previous project.

Oszi_Single_data_exchange_script: Creates a wireless communication over IP/TCP, where the Oscilloscope is sending data and afterward plotted. The communication is handled by the VISA drivers. The communication protocol, that is used for controlling the Oscilloscope, is the SCPI.

COM_PORT_REAL_TIME_with_Oszi: A serial port is open for the microcontroller to send data. The microcontroller is connected to the PC over USB. The data is then plotted in real-time. Every new package over USB is updating the plot.

Handling_COM_Port_Communication_with_Oszi: Combining both data sources (microcontroller, Oscilloscope) and plotting it the data on a single graph. Some modifications needed to be done as the sampling time and range scales have been different.

Post_Measurment_AnalyzingData: Script for comparing certain errors that have been caused due to the different sampling rates. The measurement data has been recorded at different sampling to estimate the errors.
