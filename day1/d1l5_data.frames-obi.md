---
layout: topic
title: More on data frames
author: Data Carpentry contributors
minutes: 30
---
> ## Learning Objectives
>
> *   understand the concept of a `data.frame`
> *   using sequences to index data
> *   know how to access any element of a `data.frame`


# What are data frames?

`data.frame` is the _de facto_ data structure for most tabular data and what we
use for statistics and plotting.

A `data.frame` is a collection of vectors of identical lengths. Each vector
represents a column, and each vector can be of a different data type (e.g.,
characters, integers, factors). The `str()` function is useful to inspect the
data types of the columns.


<!--- talk about colClasses argument?, row names?  --->

#### Creating ```data.frames``` manually

You can also create `data.frame` manually with the function `data.frame()`. This
function can also take the argument `stringsAsFactors`. Compare the output of
these examples, and compare the difference between when the data are being read
as `character` and when they are being as `factor`.

```
>example_data <- data.frame(
                            animal=c("dog", "cat", "sea cucumber", "sea urchin"),
                            feel=c("furry", "furry", "squishy", "spiny"),
                            weight=c(45, 8, 1.1, 0.8)
                            )
>str(example_data)
>
```
```
#> 'data.frame':    4 obs. of  3 variables:
#>  $ animal: Factor w/ 4 levels "cat","dog","sea cucumber",..: 2 1 3 4
#>  $ feel  : Factor w/ 3 levels "furry","spiny",..: 1 1 3 2
#>  $ weight: num  45 8 1.1 0.8
```

```
>example_data <- data.frame(
                            animal=c("dog", "cat", "sea cucumber", "sea urchin"),
                            feel=c("furry", "furry", "squishy", "spiny"),
                            weight=c(45, 8, 1.1, 0.8),
                            stringsAsFactors=FALSE)
                            )
>
>str(example_data)
```

```
#> 'data.frame':    4 obs. of  3 variables:
#>  $ animal: chr  "dog" "cat" "sea cucumber" "sea urchin"
#>  $ feel  : chr  "furry" "furry" "squishy" "spiny"
#>  $ weight: num  45 8 1.1 0.8
```

#### Creating ```data.frames``` by importing data with `read.table`

A `data.frame` can be created by the functions `read.table()`. In the following line we are creating the **data.frame** named "other_survey" by reading the file we just downloaded and specifying  that there is a header row with ``` header=T ``` , and that values in column are tab-delimited with ``` sep="\t" ```:

```
my_data  <- read.table ("my_file.txt", header=T , sep="\t" )
```

For more options check the [manual](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/read.table.html).


### Challenge

1. There are a few mistakes in this hand crafted `data.frame`, can you spot and
fix them? Don't hesitate to experiment!

    ```
    author_book <- data.frame(author_first=c("Charles", "Ernst", "Theodosius"),
                              author_last=c(Darwin, Mayr, Dobzhansky),
                              year=c(1942, 1970))
    ```

2. Can you predict the class for each of the columns in the following example?

    ```
    country_climate <- data.frame(country=c("Canada", "Panama", "South Africa", "Australia"),
                                   climate=c("cold", "hot", "temperate", "hot/temperate"),
                                   temperature=c(10, 30, 18, "15"),
                                   northern_hemisphere=c(TRUE, TRUE, FALSE, "FALSE"),
                                   has_kangaroo=c(FALSE, FALSE, FALSE, 1))
    ```

Check your guesses using `str(country_climate)`. Are they what you expected? Why? Why not?

R coerces (when possible) to the data type that is the least common denominator and the easiest to coerce to.


# Inspecting `data.frame` objects

We already saw how the functions `head()` and `str()` can be useful to check the
content and the structure of a `data.frame`. Here is a non-exhaustive list of
functions to get a sense of the content/structure of the data.

* Size:
    * `dim()` - returns a vector with the number of rows in the first element, and
	  the number of columns as the second element (the dimensions of the object)
    * `nrow()` - returns the number of rows
    * `ncol()` - returns the number of columns
* Content:
    * `head()` - shows the first 6 rows
    * `tail()` - shows the last 6 rows
* Names:
    * `names()` - returns the column names (synonym of `colnames()` for `data.frame`
	objects)
   * `rownames()` - returns the row names
* Summary:
   * `str()` - structure of the object and information about the class, length and
	content of  each column
   * `summary()` - summary statistics for each column

Note: most of these functions are "generic", they can be used on other types of
objects besides `data.frame`.

# Indexing and sequences

If we want to extract one or several values from a vector, we must provide one
or several indices in square brackets, just as we do in math. For instance try the following:

```
>animals <- c("mouse", "rat", "dog", "cat")
>animals[2]
>
>animals[c(3, 2)]
>animals[2:4]
>
>more_animals <- animals[c(1:3, 2:4)]
>more_animals
```

R indexes start at 1. Programming languages like Fortran, MATLAB, and R start
counting at 1, because that's what human beings typically do. Languages in the C
family (including C++, Java, Perl, and Python) count from 0 because that's
simpler for computers to do.

`:` is a special function that creates numeric vectors of integers in increasing
or decreasing order, test `1:10` and `10:1` for instance. The function `seq()`
(for sequence) can be used to create more complex patterns:

```
>seq(1, 10, by=2)
>seq(5, 10, length.out=3)
>seq(50, by=5, length.out=10)
>seq(1, 8, by=3) # sequence stops to stay below upper limit
```

Our survey data frame has rows and columns (it has 2 dimensions), if we want to
extract some specific data from it, we need to specify the "coordinates" we want
from it. Row numbers come first, followed by column numbers.

```
>surveys[1, 1]   # first element in the first column of the data frame
>surveys[1, 6]   # first element in the 6th column
>surveys[1:3, 7] # first three elements in the 7th column
>surveys[3, ]    # the 3rd element for all columns
>surveys[, 8]    # the entire 8th column
>head_surveys <- surveys[1:6, ] # surveys[1:6, ] is equivalent to head(surveys)
```

### Challenge

1. The function `nrow()` on a `data.frame` returns the number of rows. Use it,
   in conjuction with `seq()` to create a new `data.frame` called
   `surveys_by_10` that includes every 10th row of the survey data frame
   starting at row 10 (10, 20, 30, ...)

```
```

<!---
```
surveys_by_10 <- surveys[seq(10, nrow(surveys), by=10), ]
```
--->
