# RNS_LRA_EC_Scalar_Multiplier
In this repository, a Elliptic Curve scalar Multiplication algorithm Using residue Number System is implemented in C. 
The implementation is side channel attack protected using the Leakage Resistant Arithmetic (LRA) approach.
Originally, the code is evaluated in a Raspberry pi 3 embedded system but it can be generalized for any other Device.
It needs the gmp library in order to be fully functional.


The methods that are implemented in the C code library, are described in the following papers:

Fournaris, A. P., Papachristodoulou, L., & Sklavos, N. (2017, April). Secure and Efficient RNS Software Implementation for Elliptic Curve Cryptography. In Security and Privacy Workshops (EuroS&PW), 2017 IEEE European Symposium on (pp. 86-93). IEEE.
(http://ieeexplore.ieee.org/abstract/document/7966976/)

Fournaris, A. P., Papachristodoulou, L., Batina, L., & Sklavos, N. (2016, April). Residue number system as a side channel and fault injection attack countermeasure in elliptic curve cryptography. In Design and Technology of Integrated Systems in Nanoscale Era (DTIS), 2016 International Conference on (pp. 1-4). IEEE.
(https://ieeexplore.ieee.org/abstract/document/7483807/)

