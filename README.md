# MMP_2D_FT
Build instructions:
1) navigate to the Encoder directory.
2) make clean
3) make
4) navigate to the Decoder directory.
5) make clean
6) make

To encode a raw black and white 8 bits/pixel image (PGM without the header):

[path to encoder]MMP2D-FT-MP -f [image file] -lin [number of lines] -col [number of columns]
-nmax [number of lines of the largest block] -mmax [number of columns of the largest block]
-lambda [Lagrange multiplier value] -nt [number of threads to use]

The program will create a file of the same name as the input image but with appended extension .mmp

To decode an .mmp file:

[path to decoder]UNMMP2D-FT-MP -f [mmp file] -lin [number of lines] -col [number of columns]
-nmax [number of lines of the largest block] -mmax [number of columns of the largest block]
