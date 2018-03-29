# 384-bit-Rate-Pairing-Computation-in-SM9
-20180329-
-edit by LuoGuiwen-

384-bit Rate Pairing Computation in SM9

We propose a program to compute Rate Pairing in cryptography algorithm SM9. It should be used with GNU MP. You can employ this program to compute arbitrary bit length Rate pairing on BN curve in SM9 with the twist y^2=x^3+b/w^6,and the extension field describing in Chapter 4 of my paper <<Constructing 384-bit SM9 System Parmeters>>.

There are 2 C++ source files in the program,Fp12ByLGWgmp.h and ratepairingcompute3.0.cpp, and an executable file, ratepairingcompute3.0.o. 
You can directly run ratepairingcompute3.0.o under linux environment, or compile using following command,

g++ -g -O2 ratepairingcompute3.0.cpp -o ratepairingcompute3.0 -lgmpxx -lgmp -I ~/Documents/LGWProgrammingTraining/384bitratepairingcompute20180101

The command after -I is the dirname of where you store the sourcecode files.

Sincerely thanks to GNU MP library, and those academic papers I referenced in this C++ programme, which are listed in my paper.
