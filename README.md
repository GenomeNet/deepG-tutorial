# deepG Tutorials

The deepG library can be used for applying deep learning on genomic data. The library supports creating neural network architecture,
automation of data preprocessing (data generator), network training, inference and visualizing feature importance. 

For ISMB users: please jump directly to Step 1 since the VM contains demo data and preinstalled deepG.

## Overview of deepG

## 16S rRNA gene detection

### Step 0a. Install library

Before starting, make sure to have deepG installed.

```r
library(deepG)
```

### Step 0b. Download example files

This repository contains training files of 16S rRNA genes in `16s/` and a negative dataset `bacteria/`. To get this data

```bash
git clone https://github.com/GenomeNet/tutorial
```

### Setp 1. Define architecture

A _deepG_ neural network architecture consists of multiple layers. The user can specify the layer types and the size of the layers. 

In general the user can choose between:
- *LSTM layers* (long short-term memory) which are specifically designed to processes sequential data where order of data is important by using feedback connections.
- *CNN layers* (convolutional neural network) which are usually applied to images or audio data but can also be used for natural language processing or genomic sequences. Contrary to vanilla feedforward networks, they are able to process spatial relations in the data.
- The last dense layer has a softmax activation and determines how many targets we want to predict. This output gives a vector of probabilities, i.e. the sum of the vector is 1 and each entry is a probability for one class.  

To design the network we have choosen the following parameters
- `maxlen = 500` - number of nucleotides processed in one sample
- `layer_lstm = NULL` - number of LSTM cells, here only CNN are used
- `layer_dense = c(2)` - number of neurons in last layer, since it is a binary taks the output here is 2-dimensional
- `vocabulary.size = 4` - 4 letters in the DNA alphabet
- `kernel_size = c(12, 12, 12)` - size of individual CNN windows 
- `filters = c(32, 64, 64)` - number of CNN filters per layer
- `pool_size = c(3, 3, 3)` -  size of max pooling per layer

- for example the LSTM layer may be bidirectional (runs input in two ways) or stateful (considers dependencies between batches).  

The following implementation creates a model with 3 CNN layer (+ batch normalization and max pooling) followed by a LSTM layer and a final dense layer

```{r}
model <- create_model_lstm_cnn(
maxlen = 500,
layer_lstm = NULL,
layer_dense = c(2), # predict 2 classes
vocabulary.size = 4,
kernel_size = c(12, 12, 12),
filters = c(32, 64, 128),
pool_size = c(3, 3, 3),
learning.rate = 0.001
)
```
