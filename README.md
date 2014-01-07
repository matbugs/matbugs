matbugs
=======

A MATLAB interface to WinBugs.


MATBUGS is a Matlab interface for
<a href="http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/contents.shtml">WinBugs</a>
and <a href="http://www.math.helsinki.fi/openbugs/Home.html">OpenBugs</a>,
which are programs for Gibbs sampling applied to hierarchical Bayesian models.


Written by Kevin Murphy and Maryam Mahdaviani, August 2005.


*Download*
<a href="https://github.com/matbugs/matbugs/blob/master/matbugs.m">matbugs.m</a> is the only file you need.
However, there are several demos that show how to use matbugs that are worth downloading too.
Just press the 'download ZIP' button on github, or type `git clone https://github.com/matbugs/matbugs.git'.


*Usage* 

Here is a basic example of how to use matbugs:

   `[samples, stats] = matbugs(dataStruct, modelFileName, 
        		'init', initStructs,
	        	'view', 1, 'nburnin', 1000, 'nsamples', 500, 
		        'thin', 10, 
	        	'monitorParams', {'theta', 'mu_theta', 'sigma_theta'},
	        	'Bugdir', 'C:/Program Files/WinBUGS14');`

You can see that optional arguments are passed in the form 'name', value.
To use openbugs, pass in the arguments `'openbugs', 1`.

It is very important that you initialise the variables (using initStructs) to sensible
values - if BUGS samples nodes from a vague prior,  they will likely be
extreme, and the program will crash!

There are two outputs.
|samples| is a structure, containing one field per variable
that you have chosen to monitor e.g. samples.theta(c,s) is the value
of sample s in chain c. For vectors, use samples.theta(c,s,i)
and for matrices, samples.theta(c,s,i,j).

|stats| is a structure containing the mean, variance and EPSR
statistic  computed over samples for each monitored variable.


*Demos*

The examples below give more details on how to use the program

<a href="http://matbugs.googlecode.com/svn/trunk/demos/schoolsDemo/schools_writeup.html">schools</a>

<a href="http://matbugs.googlecode.com/svn/trunk/demos/biratsDemo/birats_writeup.html">birats</a>

<a href="http://matbugs.googlecode.com/svn/trunk/demos/seedsSimpleDemo/seeds_demo.m"><seeds</a>



*How it works*

Matbugs automatically generates the script file from the model (e.g.,
<a href="http://code.google.com/p/matbugs/source/browse/trunk/demos/schoolsDemo/script.txt">
see this example script file</a>, calls WinBUGS, reads in the results 
(using <a href="http://www.lce.hut.fi/research/compinf/bugsmatlab/">bugs2mat</a>),
and then performs simple convergence diagnostics. Thus it provides complete functionality, very similar to
<a href="http://cran.r-project.org/web/packages/R2WinBUGS/index.html">R2WinBUGS</a> and
<a href="http://www.stat.columbia.edu/~gelman/bugsR">bugsR</a>.

*Issues on MS Windows systems*

If you are using Windows 7 or Windows Vista, make sure you install winbugs
in a directory to which you have write permission - don't use
'C:\Program Files\Winbugs14'. I use 'C:\kmurphy\Programs\Winbugs14' and then
modify the 'Bugdir' argument when calling matbugs.


