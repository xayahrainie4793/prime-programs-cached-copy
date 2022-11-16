
PRIMO 4.3.3 LX64 for Ubuntu 18.04 64-bit

Copyright 2001-2020, Marcel Martin. 
All rights reserved. 



INSTALLING
----------

Create a directory and copy the 60 following files in this directory:

 primo
 primo.d01
 primo.d02
 primo.d03
 primo.d04
 primo.d05
 primo.d06
 primo.d07
 primo.d08
 primo.d09
 primo.d10
 primo.d11
 primo.d12
 primo.d13
 primo.d14
 primo.d15
 primo.d16
 primo.d17
 primo.d18
 primo.d19
 primo.d20
 primo.d21
 primo.d22
 primo.d23
 primo.d24
 primo.d25
 primo.d26
 primo.d27
 primo.d28
 primo.d29
 primo.d30
 primo.d31
 primo.html
 primo.x01
 primo.x02
 primo.x03
 primo.x04
 primo.x05
 primo.x06
 primo.x07
 primo.x08
 primo.x09
 primo.x10
 primo.x11
 primo.x12
 primo.x13
 primo.x14
 primo.x15
 primo.x16
 primo.x17
 primo.x18
 primo.x19
 primo.x20
 primo.x21
 readme.txt 
 stk4330
 stk4331
 stk4332
 stk4333
 verifier-f4.txt

!!! Do not forget to make executable the primo, stk4330, stk4331, stk4332, 
and stk4333 files.



RUNNING (CERTIFICATION)
-----------------------

1) Double-click on primo.
2) Open the menu, click on the "Setup..." item and select a value in the the 
   "Max number of concurrent tasks" editor (the value should be chosen 
   according to the number of cores of the processor).
3) Select the "Certification" page.
4) Select the trial division parameters and click on the "Build prime table" 
   button.
5) Load a .in file with the "Load" button.

On step #4, depending on the selected values, you may have to increase the
Linux kernel.shmmax value. With Ubuntu, you can do it with the following
command line :

 sudo sh -c 'echo SIZE > /proc/sys/kernel/shmmax'

Replace "SIZE" with the value indicated in the Primo status line.
This modification will last until Linux is rebooted (but it is possible to set
a permanent modification, see the Ubuntu doc).

