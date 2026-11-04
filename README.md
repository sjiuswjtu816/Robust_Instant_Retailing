# Data and Code for Robust Offline-Network Design and Quantity-Discount Procurement for Instant Retailing with Store Fulfillment

-This folder contains the source code, input data, and results for the numerical experiments and case study.

## Experiments and Corresponding Tables

Table(s)   Manuscript section
Table 2    6.1.1 Comparing various approaches
Table 3    6.1.2 The value of generalized budget uncertainty set
Table 4    6.2.1 The impact of multiple store classes
Table 5    6.2.2 The impact of quantity-discount procurement
Table 6    6.2.3 The impact of online-channel preference
Tables A1-A2 6.2.4 The impacts of key financial parameters

## Folder 'Data'

-The JSON datasets follow the naming convention: 'DataSet_<I>x<T>x<J>x<K>_<label>.json'.
-'Products_Prices_Data.xlsx' contains the products price data for the case study.
-The folder 'SensiAnaly-DemandDistribution' contains data for the different beta distributions considered in Section 6.1.2. In the file naming, 'leftD1, leftD2, leftD3, rightD1, rightD2, rightD3' correspond to the Beta distributions (1,9), (1,3), (2,3), (9,1), (3,1), and (3,2), respectively. If the file name contains '20-', it indicates that the support set is [20,60], otherwise, it is [30,50]. 

## Folder 'Code'

-'EVPI.py', 'DET.py', and 'RO.py' are the comparison methods in Table2, and the JSON files in the 'Data' folder are read via the 'computing.py' file.

## Folder 'Results'

-Each subfolder contains the results corresponding to the indicated table in the manuscript.
-Each file in the subfolder 'Table7' is named according to 'DRO-BD-100×4×6×6_sample-500-\vartheta_r-\vartheta_p'. 
-Each file in the subfolder 'Table8' is named according to 'DRO-BD-100×4×6×6_sample-500-\vartheta_r-\vartheta_h-\vartheta_f'.
