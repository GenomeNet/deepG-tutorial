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

### Step 1. Define the architecture

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

### Step 2. Inspect the training data and define paths

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

```
## # A tibble: 1 x 2
## Header Sequence
## <chr> <chr>
## 1 16S_rRNA::CP015410.2:15033~ TATGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAAT
```

We can output the number of files

```r
cat("number of files in 16S train:", length(list.files(path_16S_train)), "\n")

cat("number of files in bacteria train:", length(list.files(path_bacteria_train)), "\n")
```

### Step 3. Training the momdel

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

After training is finished, the evaluation will be printed to the screen

```
## Trained on 15 samples (batch_size=NULL, epochs=5)
## Final epoch (plot to see history):
## loss: 0.009258
## acc: 0.9987
## f1: 0.9987
## val_loss: 0.002073
## val_acc: 1
## val_f1: 1
## lr: 0.001
```

### Step 3. Load the model

The `trainNetwork` function saves the architecture and weights of a model after every epoch using checkpoints. The checkpoints
get stored in a binary h5 format to save space. The file names contain the corresponding epoch, loss and accuracy. For example,
we can display the checkpoints from binary classification model for 16S/bacteria.

```r
cp <- list.files(file.path(checkpoint_path, "16S_vs_bacteria_checkpoints"),
full.names = TRUE)
print(basename(cp))
```

After training, we can load a trained model and continue training or use the model for predictions.

Let’s create a model with random weights identical to our 16S/bacteria classifier and make some predictions.

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

Lets now evaluate the performance of this model on 1000 samples, 500 from each class

```r
eval_model <- evaluateFasta(fasta.path = c(path_16S_validation,
path_bacteria_validation),
model = model,
batch.size = 100,
step = 100,
label_vocabulary = c("16s", "bacteria"),
numberOfBatches = 10,
mode = "label_folder",
verbose = FALSE)

eval_model[["accuracy"]]
```

```
## [1] 0.499
```

```r
eval_model[["confusion_matrix"]]
```

```
## Truth
## Prediction 16s bacteria
## 16s 0 1
## bacteria 500 499
```

As expected, the performance is not better than random guessing. Let’s repeat evaluation but load the weights of our pretrained model

```r
weight_path <- cp[length(cp)]
model <- keras::load_model_weights_hdf5(model, weight_path)

now we run the evaluation again

```
eval_model <- evaluateFasta(fasta.path = c(path_16S_validation,
path_bacteria_validation),
model = model,
batch.size = 100,
step = 100,
label_vocabulary = c("16s", "bacteria"),
numberOfBatches = 10,
mode = "label_folder",
verbose = FALSE)
eval_model[["accuracy"]]
```

```
## [1] 0.995
```

```r
eval_model[["confusion_matrix"]]
````

```
## Truth
## Prediction 16s bacteria
## 16s 495 0
## bacteria 5 500
```

This model outperforms the random model we evaluated before. 


### Step 4. Use the model to make predictions

Once we have trained a model, we may use the model to get the activations of a certain layer and write the states to an h5 file. First, we apply our model to a file from our 16S validation set.

```r
model <- keras::load_model_hdf5(weight_path, compile = FALSE)
model <- keras::load_model_weights_hdf5(model, weight_path)
```

it is good to know the architecture designs

```r
maxlen <- model$input$shape[[2]]
num_layers <- length(model$get_config()$layers)
layer_name <- model$get_config()$layers[[num_layers]]$name
cat("get output at layer", layer_name)
``

predict a sample in the validation set

```r
fasta.path <- list.files(path_16S_validation, # make predictions for 16S file
full.names = TRUE)[1]
fasta.file <- microseq::readFasta(fasta.path)
head(fasta.file)
```

```
## # A tibble: 1 x 2
## Header Sequence
## <chr> <chr>
## 1 16S_rRNA::CP015410.2:15033~ TATGAGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCCTAAT
```

define where the file will be saved

```r
sequence <- fasta.file$Sequence[1]
filename <- file.path(dir_path, "states.h5")
```

now we do the prediction using a sliding window over the sequence with a stepsize of 1 defined with `step = 1`

```
if (!file.exists(filename)) {
    writeStates(
        model = model,
        layer_name = layer_name,
        sequence = sequence,
        round_digits = 4,
        filename = filename,
        batch.size = 10,
        step = 1)
}

and acess the file as follows

```r
states <- readRowsFromH5(h5_path = filename, complete = TRUE)
```