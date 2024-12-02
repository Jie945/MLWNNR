# MLWNNR
The code of MLWNNR and datasets. 

Dataset1, Dataset2 and Dataset3 are the three datasets used in our model.

Demo.m is our main file. You can double-click to open it and run it to obtain the relevant verification results.

The complete list of files related to this model is as follows:
![Image text](https://github.com/Jie945/MLWNNR/blob/main/Picutre/1.png)

Here is an example of a code execution diagram:
![Image text](https://github.com/Jie945/MLWNNR/blob/main/Picutre/2.png)

Here is an example of a code execution result diagram:
![Image text](https://github.com/Jie945/MLWNNR/blob/main/Picutre/3.png)

## Method of Application
The input for the model consists of three files: Disease_LncRNA_association.csv, Disease_similarity.csv, and LncRNA_similarity.csv. Disease_LncRNA_association.csv stores the known associations between LncRNA and diseases, Disease_similarity.csv stores the similarity degrees between each disease and other diseases, and LncRNA_similarity.csv stores the similarity degrees between each LncRNA and other LncRNAs. The data preparation for these three files is described in Section 2.1. Materials.

The output of the model code is a two-dimensional matrix filled with LncRNA-disease association prediction scores, from which the association degree of each LncRNA with all diseases can be queried.

If you want to make predictions based on your own data, you need to create three two-dimensional matrices, convert them into .csv files, and then replace the original Disease_LncRNA_association.csv, Disease_similarity.csv, and LncRNA_similarity.csv files. In the first two-dimensional matrix, each data point represents the association between an LncRNA and a disease, with a value of 1 indicating a confirmed association and 0 indicating no confirmed association. In the second matrix, each data point represents the similarity degree between each LncRNA and other LncRNAs, with the calculation method provided in Section 2.1.2. LncRNA expression similarity. In the third matrix, each data point represents the similarity degree between each disease and other diseases, with the calculation method detailed in Section 2.1.4. Disease semantic similarity. If you have any questions, feel free to contact me via email, and I'd be happy to help you.




