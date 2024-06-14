#QRS Classifier
This repository contains the implementation of a novel QRS classifier. The classifier uses three (or more) QRS detection algorithms running in parallel, processes the timestamps of each detection using a decision tree, and returns the classification results.

**Database**

The code was evaluated and tested using the MIT-BIH Arrhythmia Database, which is available here: https://physionet.org/content/mitdb/1.0.0/

**Installation and Setup**

1. Clone this repository:
    ```bash
    git clone [link here]
    cd novel-qrs-classifier
    ```


**Running the Code**

To run the code, execute the following Jupyter Notebooks in order:

1. **Database Preprocessing**: 
    - `1_Database_preprocessing.ipynb`
    - This notebook produces two pandas dataframes and saves them as pickle files in the `dataframes` folder.

2. **Input to Classifier**:
    - `2_Input_to_classifier.ipynb`
    - This notebook processes the previously created dataframes further and produces a third dataframe.

3. **Classifier**:
    - `3_Decision_Tree.ipynb`
    - This notebook uses the final dataframe to train and evaluate the classifier.

**Additional Resources**

The folder `Other_notebooks` contains some useful code for visualizing the database and dataframes.

**Credits**

The MIT-BIH Arrhythmia Database is an essential part of this project. More information about the database can be found here: https://physionet.org/content/mitdb/1.0.0/

**Feedback and Questions**

If you have any questions or feedback, feel free to reach out to [jeremayaa](https://github.com/jeremayaa).