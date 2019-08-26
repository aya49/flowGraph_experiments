there exists
- relationship between cells: overlap of markers; "parent/child/sibling"
- relationship between edges: 



compare with: Flowtype and cytodx

7 —> 4 pages abstract draft out

- Greatest mean separation in genetech
- my tube3: CD34+CD117+CD56-
- blasts tube4: CD13+CD16-cd11c-

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