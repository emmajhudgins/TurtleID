# Turtle ID

Scripts for use at ISTS 2019 in concert with Jillian A. Hudgins' presentation. Examples of RMark and RCapture analyses on turtle photo ID data.

Written by Emma J. Hudgins.
For questions, e-mail emma.hudgins@mail.mcgill.ca

The scripts contained here use RCapture and RMark packages. 

RCapture installation proceeds like any R package installation, while RMark requires additional installation and configuration of the MARK program.

For instructions to install MARK (non-trivial especially on Mac, 
requires XCode command line tools, Homebrew, comfort using the Terminal), see:
http://www.phidot.org/software/mark/downloads/index.html
If you are having trouble with MARK, robust design methods are available through *RCapture::robustd.t* and *RCapture::robustd.0*

For an overview of RCapture methodology, see:
https://www.jstatsoft.org/article/view/v019i05
as well as the package description on CRAN. Bear in mind that while RCapture estimates a trap effect, the trap effect model cannot output demographic information, so the estimates are derived from the model without the trap effect.

For an overview of Robust Design methodology, see:
http://www.phidot.org/software/mark/docs/book/pdf/chap15.pdf

Sample data are a subset of Maldives Turtle Photo ID records collated by Jillian A. Hudgins
