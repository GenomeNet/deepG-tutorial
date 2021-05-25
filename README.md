# deepG Tutorials

The deepG library can be used for applying deep learning on genomic data. The library supports creating neural network architecture,
automation of data preprocessing (data generator), network training, inference and visualizing feature importance. 

For ISMB users: please jump directly to Step 1 since the VM contains demo data and preinstalled deepG.

## Overview of deepG

(TODO: add image here)

## 16S rRNA gene detection

### Step 0a. Install the library

Before starting, make sure to have deepG installed.

```r
library(deepG)
```

### Step 0b. Download the example files

This repository contains training files of 16S rRNA genes in `16s/` and a negative dataset `bacteria/`. To get this data

```bash
git clone https://github.com/GenomeNet/tutorial
```

### Setp 1. Define the architecture

A _deepG_ neural network architecture consists of multiple layers. The user can specify the layer types and the size of the layers. 

In general the user can choose between:
- *LSTM layers* (long short-term memory) which are specifically designed to processes sequential data where order of data is important by using feedback connections.
- *CNN layers* (convolutional neural network) which are usually applied to images or audio data but can also be used for natural language processing or genomic sequences. Contrary to vanilla feedforward networks, they are able to process spatial relations in the data.
- A last *dense layer* has a softmax activation (Softmax converts a vector of values to a probability distribution) and determines how many targets we want to predict.

To design the network we have choosen the following parameters
- `maxlen = 500` - number of nucleotides processed in one sample
- `layer_lstm = NULL` - number of LSTM cells, here only CNN are used
- `layer_dense = c(2)` - number of neurons in last layer, since it is a binary taks the output here is 2-dimensional
- `vocabulary.size = 4` - 4 letters in the DNA alphabet
- `kernel_size = c(12, 12, 12)` - size of individual CNN windows 
- `filters = c(32, 64, 64)` - number of CNN filters per layer
- `pool_size = c(3, 3, 3)` -  size of max pooling per layer

The following implementation creates a model with 3 CNN layer (+ batch normalization and max pooling) followed by a final dense layer

```r
model <- create_model_lstm_cnn(
maxlen = 500,
layer_lstm = NULL,
layer_dense = c(2),
vocabulary.size = 4,
kernel_size = c(12, 12, 12),
filters = c(32, 64, 128),
pool_size = c(3, 3, 3),
learning.rate = 0.001)
```

### Setp 2. Inspect the training data and define paths

Input data must be a collection of files in FASTA or FASTQ format. When training a binary or multi-class model the esiest option is to put training files belonging to a class to a seperate folder. It is recommendet to split the dataset to train/validation to measure overfitting.

```r
path <- "/home/rmreches/tutorial"
path_16S_train <- file.path(path, "16s/train")
path_16S_validation <- file.path(path, "16s/validation")
path_bacteria_train <- file.path(path, "bacteria/train")
path_bacteria_validation <- file.path(path, "bacteria/validation")
```

Furthermore we need file paths where we store logs

```r
checkpoint_path <- file.path(path, "checkpoints")
tensorboard.log <- file.path(path, "tensorboard")
dir_path <- file.path(path, "outputs")
if (!dir.exists(checkpoint_path)) dir.create(checkpoint_path)
if (!dir.exists(tensorboard.log)) dir.create(tensorboard.log)
if (!dir.exists(dir_path)) dir.create(dir_path)
```

Lets take a peek into the training dataset and output the first sequence

```r
print(microseq::readFasta(list.files(path_16S_train, full.names = TRUE)[1]))
```

We can output the number of files

```r
cat("number of files in 16S train:", length(list.files(path_16S_train)), "\n")

cat("number of files in bacteria train:", length(list.files(path_bacteria_train)), "\n")
```

### Setp 3. Training the momdel

We need to specify our settings to the `trainNetwork` function, especially
- `train_type` - set to `"label_folder"` since we have a folder per class
- we supply the path to training/validation data with `path` and `path.val`
- `labelVocabulary = c("16s", "bacteria")` - the label names in same order as `path` and `path.val`
- `model = model` - here we supply the model created in *Step 1*
- `batch.size = 256` the number of samples that are processed in parallel (half of batch is 16S and other half bacteria data)
- `steps.per.epoch = 15` specifies that 15 batches will be processed before a epoch is finished
- `epochs = 5` - training will be finished after 5 epochs
- `proportion_per_file = c(0.95, 0.05)` - randomly select 95% of 16S and 5% of bacteria files, since bacteria files are much larger

```r
trainNetwork(train_type = "label_folder",
model = model,
path = c(path_16S_train, path_bacteria_train),
path.val = c(path_16S_validation path_bacteria_validation),
labelVocabulary = c("16s", "bacteria"), # label names
checkpoint_path = checkpoint_path,
tensorboard.log = tensorboard.log,
validation.split = 0.2,
run.name = "16S_vs_bacteria",
batch.size = 256,
steps.per.epoch = 15,
epochs = 5,
save_best_only = FALSE,
step = c(100, 500), # smaller step size for 16S
output = list(none = FALSE,
    checkpoints = TRUE,
    tensorboard = TRUE,
    log = FALSE,
    serialize_model = FALSE,
    full_model = FALSE),
tb_images = TRUE,
proportion_per_file = c(0.95, 0.05))
```