# bambi
A set of programs to manipulate SAM/BAM/CRAM files, using HTSLIB

to install:


clone into a local directory (eg 'bambi')

`cd bambi`

`autoreconf -i`

`configure --prefix=<place to install> --with-htslib=<directory containing a htslib release>`

(eg configure --prefix=/usr/local/bambi --with-htslib=/usr/local/htslib/1.3)

`make`

`make check` (to run tests)

