# Code for "Mapping the Extinction and Reddening of the Near Infrared Using UKIRT"

Public cleaned version of code used to create extinction maps.

This code allows one to create maps towards the Galactic Center using UKIRT photometry data (CASU or PSF) for all years. In order to run the code, make sure you check the options file and edit any file directories to match those which you will use.

Order to run code (assuming all works):
1. readData.py
2. converter.py
3. mapper.py (doesn't technically have to be ran if choiceMap.py is used)
4. make_map.py or choiceMap.py (either can be ran, choiceMap just consolidates changes)