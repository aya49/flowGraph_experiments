# flowtype_metrics

code adapted from flowtype_IMPC + _FlowCAP; _func.R adapted from 2dgate 20190329

data sets
- impc 
  - class: gene (wildtype wt control gene & various knockout ko genes), gender; files are normalized based on ko genes
- flowcap-II aml
  - class: aml (half of healthy patients are used as controls, healthy, aml)
- pregnancy 
  - class: 4 time points of pregnancy = early, mid, late, 6 weeks postpartum
- genentech
  - class: bone marrow and whole blood mixed and pure
  - 2018 results based on 3 donors (1 is an outlier) greatest mean separation:
    - my tube3: CD34+CD117+CD56-
    - blasts tube4: CD13+CD16-cd11c-
  - 2019-11-19 results based on a single random donor D
    - singlets tube4: CD13+CD11c-
    - singlets tube4 (suggested experimentally by GNE): CD13+CD16-


scripts
- dependency between scripts are indicated by their number e.g. 02 needs to be ran before 03
- scripts ending in - (working or on hold) or _ (old) are not run

methods described in: https://www.overleaf.com/7987944831tjbjkzvpmbnt

see script comments for which methods in the above document are ran

compare with: Flowtype and cytodx

Instead of overlap: number of standard deviations of control between means

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