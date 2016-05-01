---
layout: topic
title: Before we start
author: Data Carpentry contributors
minutes: 30
---


------------

> ## Learning Objectives
>
> * Articulating motivations for this lesson
> * Set up participants to have a working directory with a `data/` folder inside
> * Introduce R syntax
> * Point to relevant information on how to get help, and understand how to ask well formulated questions

------------

# Before we get started

## What is R?

R is a free language and environment for statistical computing and graphics. You can search all information about [R on the R project website](http://www.R-project.org)

R system contains two major components:
* Base system – contains the R language software and the high priority add-on packages listed on site.
* User contributed add-on package


## Install R

Go to the [download page](http://CRAN.R-project.org), choose the version that suits your computer operating system and follow the instructions.


## Launch R

If you are using a computer with Linux or Mac, open terminal and launch R by typing "R".

If you are using a computer with Windows, launch the application from the desktop icon.   

In all cases you will see something like this:

ADD IMAGE


the symbol `>` is the prompt, the place where you will communicate with your machine using R. All instructions should be written here either by typing or copy-pasting.

By pressing `Enter` after an instruction, R will try to execute it, and when ready, show the results and
come back with a new `>` prompt to wait for new commands.

Instructions are given as commands. A typical command looks like:

```
> round(3.1415926535897932384626433832795028841971)
[1] 3
>
```

If R is still waiting for you to enter more data because it isn't complete yet,
the console will show a `+` prompt. It means that you haven't finished entering
a complete command. This is because you have not 'closed' a parenthesis or
quotation. If this happen press `Ctr-c`,  this should help you out of trouble.

```
> round(3.1415926535897932384626433832795028841971
+
+
>
```

## Where am I? The workspace

To learn in which part of the machine you are use the command `getwd()`

```
>getwd()
>getwd()
[1] "/home/username/userfolder"
>
```
To see the content of the folder use the command `dir()`

```
>dir()
[1] "somefile" "someotherfile"
>
```

## Organizing your working directory

You should separate the original data (raw data) from intermediate datasets that
you may create for the need of a particular analysis. For instance, you may want
to create a `data/` directory within your working directory that stores the raw
data, and have a `data_output/` directory for intermediate datasets and a
`figure_output/` directory for the plots you will generate.

To create a directory inside R use the command `dir.create`

```
> dir.create("mydata")
> dir()
[1] "mydata"    "somefile" "someotherfile"   
>
```

If you want to move in the new created directory use the command `setwd`:

```
> getwd()
[1] "/home/username/userfolder"
>
>dir()
[1] "mydata"    "somefile" "someotherfile"
>
>setwd("/home/username/userfolder/mydata")
dir()
[1]
>
>setwd("/home/username/userfolder/")
>
>
```

##  Quit R
Type `q()` to quit


# Basics of R

R is a versatile, open source programming/scripting language that's useful both
for statistics but also data science. Inspired by the programming language S.

* Open source software under GPL.
* Superior (if not just comparable) to commercial alternatives. R has over 7,000
  user contributed packages at this time. It's widely used both in academia and
  industry.
* Available on all platforms.
* Not just for statistics, but also general purpose programming.
* For people who have experience in programmming: R is both an object-oriented
  and a so-called [functional language](http://adv-r.had.co.nz/Functional-programming.html)
* Large and growing community of peers.

## Commenting

Use `#` signs to comment. Comment liberally in your R scripts. Anything to the
right of a `#` is ignored by R.

## Assignment operator

`<-` is the assignment operator. It assigns values on the right to objects on
the left. So, after executing `x <- 3`, the value of `x` is `3`. The arrow can
be read as 3 **goes into** `x`.  

You can also use `=` for assignments but not in
all contexts so it is good practice to use `<-` for assignments. `=` should only
be used to specify the values of arguments in functions, see below.

An example of assignment is:

`a<-4`


## Functions and their arguments

Functions are "canned scripts" that automate something complicated or convenient
or both.  Many functions are predefined, or become available when using the
function `library()` (more on that later). A function usually gets one or more
inputs called *arguments* and is written in the form of:

function (arguments)

Functions often (but not always) return a *value*. A
typical example would be the function `sqrt()`. The input (the argument) must be
a number, and the return value or output is the square root of that
number.

Executing a function ('running it') is called *calling* the function. An
example of a function call is:

`b <- sqrt(a)`

Here, the value of `a` is given to the `sqrt()` function, the `sqrt()` function
calculates the square root, and returns the value which is then assigned to
variable `b`. This function is very simple, because it takes just one
argument.

The return 'value' of a function need not be numerical (like that of `sqrt()`),
and it also does not need to be a single item: it can be a set of things, or
even a data set. We'll see that when we read data files in to R.

**Arguments** can be anything, not only numbers or filenames, but also other
objects. Exactly what each argument means differs per function, and must be
looked up in the documentation (see below). If an argument alters the way the
function operates, such as whether to ignore 'bad values', such an argument is
sometimes called an *option*.

Most functions can take several arguments, but many have so-called **defaults**.
If you don't specify such an argument when calling the function, the function
itself will fall back on using the *default*. This is a standard value that the
author of the function specified as being "good enough in standard cases". An
example would be what symbol to use in a plot. However, if you want something
specific, simply change the argument yourself with a value of your choice.

Let's try a function that can take multiple arguments `round`.

```
>round(3.14159)
[1] 3

```

We can see that we get `3`. That's because the default is to round
to the nearest whole number. If we want more digits we can see
how to do that by getting information about the `round` function.
We can use another function `args()`:

```
>args(round)
> args(round)
function (x, digits = 0)
NULL
>
```

We can also look at the
help for this function using `?round`.
```
>?round
>

```
to **quit** the manual page type `q`


We see that if we want a different number of digits, we can
type `digits=2` or however many we want.

```
>round(3.14159, digits=2)
```

If you provide the arguments in the exact same order as they are defined you don't have to name them:

```
>round(3.14159, digits=2)
```

However, it's usually not recommended practice because it's a lot of remembering
to do, and if you share your code with others that includes less known functions
it makes your code difficult to read. (It's however OK to not include the names
of the arguments for basic functions like `mean`, `min`, etc...)

Another advantage of naming arguments, is that the order doesn't matter.
This is useful when there start to be more arguments.



## Seeking help

### I know the name of the function I want to use, but I'm not sure how to use it

If you need help with a specific function, let's say `barplot()`, you can type:

```
?barplot
```

If you just need to remind yourself of the names of the arguments, you can use:

```
args(lm)
```

If the function is part of a package that is installed on your computer but
don't remember which one, you can type:

```
??geom_point
```

### I want to use a function that does X, there must be a function for it but I don't know which one...

If you are looking for a function to do a particular task, you can use
`help.search()` (but only looks through the installed packages):

```
help.search("kruskal")
```

If you can't find what you are looking for, you can use the
[rdocumention.org](http://www.rdocumentation.org) website that search through
the help files across all packages available.

### I am stuck... I get an error message that I don't understand

Start by googling the error message. However, this doesn't always work very well
because often, package developers rely on the error catching provided by R. You
end up with general error messages that might not be very helpful to diagnose a
problem (e.g. "subscript out of bounds").

However, you should check stackoverflow. Search using the `[r]` tag. Most
questions have already been answered, but the challenge is to use the right
words in the search to find the answers:
[http://stackoverflow.com/questions/tagged/r](http://stackoverflow.com/questions/tagged/r)

The [Introduction to R](http://cran.r-project.org/doc/manuals/R-intro.pdf) can
also be dense for people with little programming experience but it is a good
place to understand the underpinnings of the R language.

The [R FAQ](http://cran.r-project.org/doc/FAQ/R-FAQ.html) is dense and technical
but it is full of useful information.

### Asking for help

The key to get help from someone is for them to grasp your problem rapidly. You
should make it as easy as possible to pinpoint where the issue might be.

Try to use the correct words to describe your problem. For instance, a package
is not the same thing as a library. Most people will understand what you meant,
but others have really strong feelings about the difference in meaning. The key
point is that it can make things confusing for people trying to help you. Be as
precise as possible when describing your problem

If possible, try to reduce what doesn't work to a simple reproducible
example. If you can reproduce the problem using a very small `data.frame`
instead of your 50,000 rows and 10,000 columns one, provide the small one with
the description of your problem. When appropriate, try to generalize what you
are doing so even people who are not in your field can understand the question.

To share an object with someone else, if it's relatively small, you can use the
function `dput()`. It will output R code that can be used to recreate the exact same
object as the one in memory:

```
dput(head(iris)) # iris is an example data.frame that comes with R
```

If the object is larger, provide either the raw file (i.e., your CSV file) with
your script up to the point of the error (and after removing everything that is
not relevant to your issue). Alternatively, in particular if your questions is
not related to a `data.frame`, you can save any R object to a file:

```
saveRDS(iris, file="/tmp/iris.rds")
```

The content of this file is however not human readable and cannot be posted
directly on stackoverflow. It can however be sent to someone by email who can read
it with this command:

```
some_data <- readRDS(file="~/Downloads/iris.rds")
```

Last, but certainly not least, **always include the output of `sessionInfo()`**
as it provides critical information about your platform, the versions of R and
the packages that you are using, and other information that can be very helpful
to understand your problem.

```
sessionInfo()
```

### Where to ask for help?

* Your friendly colleagues: if you know someone with more experience than you,
  they might be able and willing to help you.
* Stackoverlow: if your question hasn't been answered before and is well
  crafted, chances are you will get an answer in less than 5 min.
* The [R-help](https://stat.ethz.ch/mailman/listinfo/r-help): it is read by a
  lot of people (including most of the R core team), a lot of people post to it,
  but the tone can be pretty dry, and it is not always very welcoming to new
  users. If your question is valid, you are likely to get an answer very fast
  but don't expect that it will come with smiley faces. Also, here more than
  everywhere else, be sure to use correct vocabulary (otherwise you might get an
  answer pointing to the misuse of your words rather than answering your
  question). You will also have more success if your question is about a base
  function rather than a specific package.
* If your question is about a specific package, see if there is a mailing list
  for it. Usually it's included in the DESCRIPTION file of the package that can
  be accessed using `packageDescription("name-of-package")`. You may also want
  to try to email the author of the package directly.
* There are also some topic-specific mailing lists (GIS, phylogenetics, etc...),
  the complete list is [here](http://www.r-project.org/mail.html).

### More resources

* The [Posting Guide](http://www.r-project.org/posting-guide.html) for the R
  mailing lists.
* [How to ask for R help](http://blog.revolutionanalytics.com/2014/01/how-to-ask-for-r-help.html)
  useful guidelines
