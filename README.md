<!-- ⚠️ This README has been generated from the file(s) "blueprint.md" ⚠️-->
[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#torchic)

# ➤ torchic

**Tools for Optimized Research in Computational HIgh-energy physiCs**

<div style="text-align: center;">
  <img src="images/torchic-logo.webp" alt="Torchic Logo" width="300" />
</div>

`torchic` is a Python package designed to support researchers in computational high-energy physics with optimized tools and utilities.

## ➤ Installation

To install `torchic`, use pip:

```bash
pip install torchic
```

## ➤ Features

Provides optimized tools specifically tailored for computational high-energy physics.
Includes a variety of modules and functions to streamline research workflows.
Accessible examples and manual tests to help validate functionality.

### `histogram` Module

The `histogram` module in `torchic` offers utilities for creating and managing histograms, including support for 1D and 2D histograms with customizable axes and loading options.

#### `AxisSpec` Class

The `AxisSpec` class provides a flexible way to define histogram axes, with parameters for binning and axis labels.

- **Attributes**
  - `bins` (int): The number of bins for the axis.
  - `xmin` (double)
  - `xmax` (double)
  - `name` (str): Name of the histogram
  - `title` (str): Title of the histogram

#### `HistLoadInfo` Class

The `HistLoadInfo` class encapsulates metadata and loading instructions for histograms, ensuring consistency in how histograms are managed and used within datasets.

- **Attributes**
  - `file_path` (str): Path to the file containing histogram data.
  - `hist_name` (str): The name of the histogram within the file.

#### `build_TH1` Function

Generates a 1D histogram (TH1) with the specified data and axis specifications.

- **Parameters**
  - `data` (array-like): The data to be binned in the histogram.
  - `axis_spec` (`AxisSpec`): The specification for the histogram's axis, including bin count, range, and label.

#### `build_TH2` Function

Generates a 2D histogram (TH2) with data mapped over two axes, suitable for representing joint distributions or correlations.

- **Parameters**
  - `data_x` (array-like): The data for the x-axis.
  - `data_y` (array-like): The data for the y-axis.
  - `x_axis_spec` (`AxisSpec`): Specification for the x-axis, including bins, range, and label.
  - `y_axis_spec` (`AxisSpec`): Specification for the y-axis, including bins, range, and label.

#### `load_hist` Function

Loads a histogram based on the specified `HistLoadInfo` metadata, providing an easy interface for accessing pre-saved histograms.

- **Parameters**
  - `load_info` (`HistLoadInfo`): An instance of `HistLoadInfo` that contains the necessary file path and options for loading the histogram.
  
The `histogram` module simplifies the process of creating, configuring, and loading histograms, supporting both 1D and 2D data analysis.


### ➤ `Dataset` Class

The `Dataset` class in `torchic` provides a structured way to handle, manipulate, and analyze datasets, especially for computational high-energy physics research. Key features include subset management, efficient access to data columns, and histogram-building capabilities.

#### `add_subset(subset_name, condition)`
Adds a new subset to the dataset. This is useful for organizing data into manageable parts or categories.

```python
dataset = Dataset()
dataset.add_subset("subset_name", condition)
```
#### `__getitem__()`
Accesses data columns. With `__getitem__`, you can retrieve either a specific column across all data or a particular column from a subset.

#### `build_hist(column_name, axis_spec, **kwargs) -> TH1F`
#### `build_hist(column_name_x, column_name_y, axis_spec_x, axis_spec_y, **kwargs) -> TH2F`
Creates a ROOT histogram from a column of the dataset.

Examples:

```python
dataset.add_subset("column_name", dataset._data["A"] < 5)
subset_column = dataset["subset_name:column_name"]
# Getting a Column from a Specific Subset

axis_spec = AxisSpec(nbins, xmin, xmax, hist_name, hist_title)
dataset.build_hist(column_name, axis_spec, subset="subset_name")
# Generates a histogram for the specified column. This is especially helpful for visualizing distributions in your data.
```

The Dataset class makes it straightforward to organize data, retrieve specific columns, and generate histograms, making it ideal for physics research workflows.


## ➤ Testing
Manual tests are provided in the tests folder to verify the functionality of each module. Before using torchic in a production environment, it’s recommended to run these tests to ensure everything works as expected.

## ➤ Examples

Refer to the examples folder for use cases and example scripts that illustrate how to integrate torchic tools into your projects.
