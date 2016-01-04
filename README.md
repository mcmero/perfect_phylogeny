# README #

### What is this repository for? ###

This code is based on the blog posts explaining the [perfect phylogeny algorithm](https://genomejigsaw.wordpress.com/2015/09/09/building-phylogenetic-trees-with-binary-traits/) and the [incomplete phylogeny algorithm](https://genomejigsaw.wordpress.com/2015/11/23/incomplete_phylogeny/), based on the corresponding papers by [Gusfield]((http://onlinelibrary.wiley.com/doi/10.1002/net.3230210104/abstract) and [Pe'er et al](http://epubs.siam.org/doi/abs/10.1137/S0097539702406510). The code takes a feature matrix, infers the presence or absence of incomplete fields, and determines whether a perfect phylogeny can be built, if so, it will output a plot of the corresponding tree. 

### How do I get set up? ###

Make sure you have these dependencies installed:

* [Numpy](http://www.numpy.org/)
* [Graphviz](http://www.graphviz.org/) (for plotting)

The algorithm can be run on a tab-separated matrix with the following format:

    .	c1	c2	c3	c4	c5
    s1	1	1	0	0	-1
    s2	0	-1	1	0	-1
    s3	1	1	0	0	1
    s4	0	0	1	1	-1
    s5	0	1	0	0	-1

The columns contain the feature names, and the rows the sample names. 1s indicate the presence of the feature in the sample, 0s indicate the absence, and -1s are indeterminate (unknown). Ensure the first column name has a dot ".", as otherwise the columns will be shifted. 

    usage: run.py [-h] [--plot] m_file [outname]

Example:

    python run.py test_matrix.txt test_tree --plot

## How do I run the unit tests? ##

    python test.py

