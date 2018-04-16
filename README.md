# About 
This is a repository for the Matlab code neceassary to conduct all analyses in our paper:

Fear Generalization as Threat Prediction: Adaptive Changes in Facial Exploration Strategies revealed by Fixation-Pattern
Similarity Analysis. [Link to bioarxiv](https://www.biorxiv.org/content/early/2017/04/15/125682)

In this paper, we used an eye-tracker to record eye-movements of humans during viewing of faces. 
We focus on how eye-movement patterns change when faces are associated with a aversive outcome.

# Initial Setup

You can download the data and the associated repositories, including this one from the Open Science Framework [here](https://osf.io/zud6h/files/).
Add these repositories to your Matlab path.
```matlab
addpath('/home/onat/Documents/Code/Matlab/FPSA_FearGen/');
addpath('/home/onat/Documents/Code/Matlab/globalfunctions//');
addpath('/home/onat/Desktop/fancycarp');
```
# Examples

###  Get the list of participants.

```matlab
>> FPSA_FearGen('get_subjects')
ans =
  Columns 1 through 45
     1     2     3     4     5     6     7     8     9    10    11    12    14    15    16    17    18    19    20    21    22    23    24    25    26    27    28    29    30    31    32    33    35    36    37    39    40    41    42    43    44    45    46    47    48
  Columns 46 through 74
    49    50    51    53    54    56    57    59    60    61    62    63    64    65    66    67    68    69    70    71    72    73    74    77    78    79    80    81    82
>> 
```

### Sanity check 1. 

In eye tracking it is usual that few trials are excluded due to blinks or bad calibration.
The following code plots the number of trials per condition across the participant pool.
The figure that is produced ensures that this is not overly unbalanced.

```matlab
FPSA_FearGen('get_trialcount',4)
```
`4` above refers to the generalization phase. (2: baseline; 3: conditioning; 4: test phase).

<img src="https://github.com/selimonat/FPSA_FearGen/blob/master/figures/01.png" height="400">

### Get the data to workspace using the `get_fixmat` action. 

Fixmat is a compact way of storing large amounts of eye movement recordings in the form of fixation points.

```matlab
>> fixmat = FPSA_FearGen('get_fixmat')
ans = 
  Fixmat with properties:

     subject: [1×118188 uint32]
       phase: [1×118188 int32]
       start: [1×118188 int32]
        stop: [1×118188 int32]
           x: [1×118188 single]
           y: [1×118188 single]
         eye: [1×118188 uint8]
    deltacsp: [1×118188 int32]
        file: [1×118188 int32]
     oddball: [1×118188 int32]
     trialid: [1×118188 int32]
         ucs: [1×118188 int32]
         fix: [1×118188 int32]
       chain: [1×118188 double]
       isref: [1×118188 double]
```
It stores every fixation's attributes in the form of separate vectors. For example, information such as x and y coordinates, participant's index, the image, the condition are stored in the Fixmat.

However, Fixmat variable above is not a simple Matlab structure though, but rather an instance of a Fixmat object, as coded in the FancyCarp Toolbox. 
The benefit of having a Fixmat object is that there are useful functions built-in in the Fixmat object, such as for example visualizing fixation density maps as heatmaps.

For example, the following code can be used to plot a fixation density map (FDM) based on all subjects (stored in `S`) during the generalization phase:
```matlab
S = FPSA_FearGen('get_subjects');
fixmat = FPSA_FearGen('get_fixmat');
v{1} = {'subject' S 'phase' 4};
fixmat.getmaps(v{:});
fixmat.plot
```
<img src="https://github.com/selimonat/FPSA_FearGen/blob/master/figures/02.png" height="300">

The `Fixmat.getmaps` method creates an FDM based on the cell array argument `v`. In the example above,`v` is used to select all fixations that belong to both phase `4` and subjects `S`. 
The method `Fixmat.plot` plots the computed FDM. In order to create a separate FDM for different conditions or participants we would create a separate cell array for each of the required filtering conditions.

```matlab
v=[];
c=0;
for ns = S([11 3 8 23]);
  c=c+1;
  v{c} = {'subject' ns 'deltacsp' 0 'phase' 4};
end
fixmat.getmaps(v{:});
fixmat.plot
```
This will plot a separate FDM for the 4 different participants. Note how different participants scan faces differently. This is the basis for many recent reports that analyzed the scan-path idiosyncrasy during viewing of faces.
And the major reason for us to come up with the FPSA as a methodology to investigate scanning eye-movement strategies humans deploy during viewing of faces.

<img src="https://github.com/selimonat/FPSA_FearGen/blob/master/figures/03.png" height="300">

### Get FDM for the Fixation-Pattern Similarity Analysis

FPSA analysis is conducted on single-participants at a time. The `get_fixmap` action can be used to gather the required FDMs. 
It is basically a wrapper around the `Fixmat.getmaps` method.

```matlab
>> maps = FPSA_FearGen('get_fixmap',fixmat,{'subject' 2});;
>> size(maps)
ans =
      250000          16
```
As you can see in `maps` each FDM is represented as a column vector. As we have 2 phases and 8 different conditions (faces), this amounts to 16 different vectors.
This representation is appropriate for the FDM similarity analysis.

### Fixation-Pattern Similarity Analysis 

The following command will run a correlation-based similarity analysis on the FDMs of single participants, before (phase = 2) and after (phase =4) aversive learning.
`{'fix',1:100}` indicates that the analysis will include 1st to 100th fixations, that is all the fixations.

`1:3` ensures that similarity analysis is ran separately for the 3 runs of the generalization phase (phase =4). This is important because the baseline phase (phase = 2) entails only one single run. In order to have a valid comparison of the similarity values across the baseline and test phases, the FDMs should contain on average same number of fixations. Computing separate FPSA matrices for each run ensures this.

```matlab
sim = FPSA_FearGen('get_fpsa_fair',{'fix',1:100},1:3);
>> sim
sim = 
  struct with fields:

    correlation: [74×120 double]
```

<img src="https://github.com/selimonat/FPSA_FearGen/blob/master/figures/04.png" height="300">

The resulting matrix `sim.correlation` contains in each row the similarity data from a given participant. 
The pair-wise similarity values between the 16 FDMs are stored across the columns in a non-redundant format.
You can visualize the resulting 16x16 matrix:

```matlab
imagesc(squareform(mean(sim.correlation)),[0.3 1.2]);
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
axis square;
colorbar;
title('Average Dissimilarity Matrix')
```

<img src="https://github.com/selimonat/FPSA_FearGen/blob/master/figures/05.png" height="300">

The first quadrant of this matrix shows the similarity relationships between the 8 FDM recorded during the baseline period, before aversive learning has taken place.
The second quadrant shows the same after the learning.










