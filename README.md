# AVC_T62U2 
*Code used for track analysis in the EURADOS uncertainty exercise*

**Revision history**:

### *********** Version V6 ******************************************
**23-JUL-2025** Fixed bug in this context with the last track. 

**17-JUN-2025** Added output of track summaries. 
### *********** Version V5 ******************************************
**16-JUN-2025** Corrected several bug in the scoring of the totals and the output (introduced on 10-MAR-2025). 

**11-JUN-2025** Introduced workaround for segmentation fault caused by the high number of events in PHITS data at low energy. 

**10-JUN-2025** Corrected a bug with the newly introduced arrays of energy imparted, total number of ions and clusters track. Currently not yet output.

**10-MAR-2025** Changed scoring of clusters per track and globally to up to the maximum ICS to enable output of more detailed ICSDs.

### *********** Version V4 *******************************************
**09-MAR-2025** Added headers for version number in version history. Reduced output to the last radial distance with entries.

**07-MAR-2025** Fixed handling data files that contain only ionizations, in which the number of saved events deviates from the number of 
simulated events plus the case of Geant4-DNA Opt4 15 eV which contain data from two simulation runs with the same event labelling. 

**06-MAR-2025** Added fix for handling of PTra data files that deviate from expected structure: They have energy loss in 5th column and energy deposit in a 7th column. 

**20-FEB-2025** Bug fix: labelling of columns in the header of output file was incorrect. 

### *********** Version V4 *******************************************
**09-MAR-2025** Added headers for version number in version history. Reduced output to the last radial distance with entries.

**07-MAR-2025** Fixed handling data files that contain only ionizations, in which the number of saved events deviates from the number of 
simulated events plus the case of Geant4-DNA Opt4 15 eV which contain data from two simulation runs with the same event labelling. 

**06-MAR-2025** Added fix for handling of PTra data files that deviate from expected structure: They have energy loss in 5th column and energy deposit in a 7th column. 

**20-FEB-2025** Bug fix: labelling of columns in the header of output file was incorrect. 

### *********** Version V3 *******************************************
**16-MAY-2024** Bug fix: When moving the first point of a new track to the top of the list after processing the previous track, the flag for ionization was omitted. (Bug reported by Pavel.)

### *********** Version V2 *******************************************
**28-JAN-2024** As decided during the last meeting of the task group on 22-Jan-2024, the code has been modified to read an additional flag indicating impact ionization (1) or electronic excitation with or without subsequent auto-ionization (2). At present this is just a provision for the final analysis to take place after a common sharepoint will be available (again). For the distributed analysis by participant, ICSD are determined scoring only impact ionizations.

### *********** Version V1 *******************************************
**07-DEC-2023** Bug fix: When moving the first point of a new track to the top of the list after processing the previous track, the value of E was not moved along. (Bug reported by Zine.)

**11-NOV-2023** Fixed several bugs based on problems reported by Zine.
1. TREIMP and TRION were not normalized to NITER (in DO_AVC)
2. NPTS was not properly set for the last track.
3. In subroutine SORT62, the case that the input data are already ascending had not been properly handled. 
4. E was not read for the first point in a track.
5. Radial scoring of the target center was flawed (the variable RSITE2 was used instead of RSITE, i.e. the square of the site radius instead of the radius.)

Further changes: 

- parameter NPMAX has been increased to 10000. 
- NP has been renamed to NPTS.
- In the main loop, the variable NPTS is now used 
  directly instead of I. 

**22-OCT-2023** Bug fix: see inline comments starting with the revision date in the corresponding code section below, at the end of the 'Process_this_track' loop. 

**04-OCT-2023** First released version
