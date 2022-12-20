## R CMD check results

0 errors | 0 warnings | 2 notes

Comments on "checking R code for possible problems" Note:
This note is triggered by the code we have accessing an external
resource known as scType. Currently, the scType tool can only be
accessed via the "source" command, so we couldn't find a good
workaround for this method of importing the tool.

Comments on "checking CRAN incoming feasibility" Note:
Our package requires the use of two packages from non-standard repos.
Therse are scPred and SingleR. A vital part of the tools functionality,
we provide information on how to install these packages on our github,
but if they cannot be installed, our tool can run without them.

* This is a new release.
