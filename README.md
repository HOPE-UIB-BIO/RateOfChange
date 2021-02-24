# RateOfChange
Script for manuscript "Rate-of-change analysis in palaeoecology revisited: a new approach" 

Abstract:
Dynamics in the rate of compositional change beyond the time of human observation are uniquely preserved in palaeoecological sequences of terrestrial and marine origin. Periods with a high rate of change (RoC) have been assigned to reflect changes due to external factors such as climate and human forcing. However, changes in sedimentation rates and sampling strategies can result in an uneven distribution of time intervals within stratigraphic sequence, which are known to affect RoC estimates and thus derived conclusions. Several approaches have been used to account for these issues (e.g. binning of levels, data smoothing, etc.) but despite the importance and frequency of these, there has been relatively little exploration of  the implications of these challenges into quantifying RoC in palaeoecology.
Here, we introduce R-Ratepol – an easy-to-use R package – that provides a new, robust numerical technique for detecting and summarising RoC patterns in complex multivariate time-ordered stratigraphical sequences. First, we compare the performance of common methods of estimating RoC and detection of periods of high RoC (peak-point) using simulated pollen-stratigraphical data with known patterns of compositional change and temporal resolution. Then, we apply our new methodology to four representative European pollen sequences to detect potential drivers of compositional change. 
We show that the successful detection of known shifts in compositional changes identified by RoC depends heavily on the smoothing methods and dissimilarity coefficients used, and the density of stratigraphical levels. Identifying the need for an improved numerical approach, we propose a new method of binning with a moving window in combination with a generalised additive model for peak-point detection. Our method shows a 68% increase in the correct detection of peak-points compared to the more traditional way of peak selection by individual levels, as well as achieving a reasonable compromise between Type I and Type II errors. Likewise, by using the four pollen sequences from Europe, we show that our approach also performs well in detecting periods of significant compositional change during known onsets of human activity, early land-use transformation, and changes in fire frequency. 
Expanding the approach using R-Ratepol to the increasing number of open access paleoecological datasets in global databases, such as Neotoma, will allow future palaeoecological and macroecological studies to quantify major changes in biotic composition or in sets of abiotic variables across broad spatio-temporal scale. 

Reference: Ondřej Mottl, John-Arvid Grytnes, Alistair W.R. Seddon, Manuel J. Steinbauer, Kuber P. Bhatta, Vivian A. Felde, Suzette G.A. Flantua, H. John B. Birks. Rate-of-change analysis in palaeoecology revisited: a new approach bioRxiv 2020.12.16.422943; doi: https://doi.org/10.1101/2020.12.16.422943
