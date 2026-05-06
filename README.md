Cost-free lifespan extension through reduced adulthood insulin-signalling is robust to natural environmental challenges
--------------------------------------------------------------------------------------------------------------------------------
Sara D. Irish*, Annabel Kimberley, Simone Immler, Alexei A. Maklakov <br>
  University of East Anglia, Norwich, UK <br>
  *Corresponding author: s.irish@uea.ac.uk <br>
---------------------------------------------------------------------------------------------------------------------------------
<br>
All data and code included here for manuscript entitled: ‘Cost-free lifespan extension through reduced adulthood insulin-signalling is robust to natural environmental challenges.’ Datasets can be found in the ‘data’ folder and all code for plots and analysis can be found in the ‘code’ folder. All data were collected at University of East Anglia, Norwich, UK. The purpose of this project was to test the whether the cost-free lifespan benefits reduced insulin signalling (IIS) in adulthood in Caenorhabditis elegans nematodes was robust to natural temperature cycles, rather than standard laboratory conditions. We used RNAi to knockdown the daf-2 gene, the IIS receptor in C. elegans (N2) after sexual maturation and exposed nematodes to daily temperature cycles between 10 and 15C, similar to what would be experienced on an autumn day in their natural habitat in Bristol, UK. We assessed their lifespan and daily reproductive output in this environment. Additionally, we performed these experiments in skn-1 and daf-16 loss of function mutants to try to find the mechanisms behind lifespan extension in semi-natural temperature conditions. 
<br>

---------------------------------------------------------------------------------------------------------------------------------

File descriptions:<br>
<br>
File: ‘Bristol_LS_mut.csv’<br>
This file contains lifespan data for C. elegans N2 wildtype strain and skn-1 and daf-16 loss of function mutant strains.<br>
Variables:<br>
Plate.ID: Nematodes were housed in plates of 10 for lifespan experiments to prevent overcrowding and starvation during the reproductive period. Plate.ID is an ID number for the plate each nematode was housed on.<br>
Worm: An individual identification number (1-10) for each worm per plate, given when the individual died. <br>
Strain: Represents the strain (N2, skn-1 (SKN), or daf-16 (DAF16)) of the individual<br>
Treatment: Represents the treatment given to each individual (EV = control provided with an empty plasmid or DAF-2 = daf-2 RNAi treatment)<br>
Death.Date: Date of death, give as DD-Month<br>
Cause: Letter representing the cause of death (D = natural death, M = matricide (when eggs hatch inside the nematode), E = explosion (when nematode cuticle fails and intestines burst through), W = walled (individual climbs onto side of plate and desiccates), or L = lost due to human error. <br>
D1: Date the experiment began<br>
Age: Number of days the individual survived (Death.Date – D1)<br>
Event: Binary variable indicating whether the individual died naturally (1, when cause = D or M) or otherwise (0, when cause = W or L)<br>
Exp: Binary variable indicating whether an individual exploded (1) or not (0)<br>
---------------------------------------------------------------------------------------------------------------------------------<br>

File: ‘Bristol_repro_mut.csv’<br>
This file contains daily reproductive output data for C. elegans N2 wildtype strain and skn-1 and daf-16 loss of function mutant strains.<br>
Variabes: <br>
Strain: Represents the strain (N2, skn-1 (SKN), or daf-16 (DAF16)) of the individual<br>
Treatment: Represents the treatment given to each individual (ev = control provided with an empty plasmid or daf = daf-2 RNAi treatment)<br>
ID: Individual ID number<br>
D1-D15: Column for each day (Days 1 through 15) of reproduction, given as a number of offspring produced each day. If NA, the individual was lost or died before the final day of reproduction. If the individual experienced death due to matricide, NA’s were later changed to 0<br>
Matricide: Binary variable indicating whether an individual died due to matricide (1) or not (0). <br>
---------------------------------------------------------------------------------------------------------------------------------<br>

File: ‘Bristol_temp.csv’<br>
This file contains the hourly temperature fluctuations over each 24-hour period of the experiment.<br>
Variables: <br>
Time: hour of day (00:00 – 23:00)<br>
Temp: temperature (in degrees Celsius) at any given hour

