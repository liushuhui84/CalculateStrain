## Strain Calculation Implementation is divided into two steps:

### 1.Atomic Column Recognition

We use Python language to train a deep learning network (U-Net) and perform batch segmentation of the raw images. The output includes the segmentation results and the corresponding atomic coordinates in an Excel file. Please refer to the `Atom-Identification.ipynb` file for the entire process.

### 2.Strain Calculation

Using the atomic column coordinates output from the first step, we perform the strain calculation. We used MATLAB to implement this, and the specific code can be found in the `strain.m` file.
