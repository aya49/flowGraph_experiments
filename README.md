# flowtype_metrics

## scripts

run/see 00_main.Rmd script(s)
- scriptsadapted from flowtype_IMPC + _FlowCAP; _func.R adapted from 2dgate 20190329
- dependency between scripts are indicated by their number e.g. 02 needs to be ran before 03
- scripts ending in - (working or on hold) or _ (old) are not run
- includes comparison with flowType count/prop featres

see 01_cytodx.R for comparisons with cytoDX

## data sets

- impc (ON-HOLD)
  - class: gene (wildtype wt control gene & various knockout ko genes), gender; 
samples should be normalized based on ko genes

- flowcap-II aml [@aghaeepour2013critical]
  - class: 43 aml (acute myeloid leukemia), 316 healthy subjects bone marrow or blood x 7 panels 
(1st panel is a control).
  - we use the 6th panel with markers $SS$, $FS$, $HLA-DR$, $CD117$, $CD45$, $CD34$, and $CD38$.
  - each sample has ~ 60,000 cells.
  - $CD34+$ increases in aml subjects.

- genentech
  - class: bone marrow and whole blood mixed and pure
  - 2019-11-19 results based on a single random donor D
    - singlets tube4: CD13+CD11c-
    - singlets tube4 (suggested experimentally by GNE): CD13+CD16-
  - earlier results based on 3 donors (1 is an outlier) greatest mean separation:
    - my tube3: CD34+CD117+CD56-
    - blasts tube4: CD13+CD16-cd11c-

- pregnancy [@aghaeepour2017immune]
  - class: 4 time points of pregnancy, early, mid, late, 6 weeks postpartum x 18 and 10 women of the training and validation cohort
  - analyzed on 13 markers ($CD123$, $CD14$, $CD16$, $CD3$, $CD4$, $CD45$, $CD45RA$, $CD56$, $CD66$, $CD7$, $CD8$, $Tbet$, $TCRgd$)
  - each sample has ~ 300,000 cells
  - Since discrepancies between subjects is a major batch effect, we further normalize for subject. For each subject, we take the calculated feature values of all of her samples and extract the difference between them and their mean.

- bodenmiller CyTOF [@bodenmiller2012multiplexed]
  - human peripheral blood from 8 x 2 BCR-XL 
(b cell receptor / Fc receptor cross linker) un/stimulated healthy subjects
  - 10 markers


## etc

NOTE: 20200225 `devtools::run_summary()` doesn't work on manual `fg_get_summary_index`.

A summary metric suggestion: number of standard deviations of control between means

apriori
- confidence: if A then B = conf(A>B) = AB/A 
- lift: if A then B = (A|B)/AB
  - ratio of the observed support to that expected if X and Y were independent
  - 1 = independant
  - >1 = dependant
  - <1 = substitues; presence of one item has negative effect on presence of other item and vice versa
- conviction: if A then B = (1-B)/(1-conf(A>B))
  - the ratio of the expected frequency that X occurs without Y (that is to say, the frequency that the rule makes an incorrect prediction) if X and Y were independent divided by the observed frequency of incorrect predictions
  - e.g. 1.2 = incorrect 20% more often (1.2 times as often) if the association between X and Y was purely random chance
- maximal: no children are >frequent
- closed: all children have smaller count than it