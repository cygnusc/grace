*** Software requirements ***
Windows 7 or later OS with DirectX 11 or newer support. If your computer is purchased on or after 2010, it probably can run this software.

Gnuplot for windows(http://sourceforge.net/projects/gnuplot/files/gnuplot/) must be installed to visualize running results. Download the latest version (gp466-win32-setup.exe or newer) for best results.

*** How to compile from the source code ***
The pre-compiled binary code has been provided for the user's convenience. To compile the code on your own, Windows 7 or later OS with Visual Studio 2012 or later are needed. 

Download the souce code from https://github.com/cygnusc/grace/tree/master/Grace.0.01.source and double click the magnetic4.sln to open the visual studio C++ project. 

You will also need C++ AMP fft library from https://ampfft.codeplex.com/releases/view/92384. Put the .lib, .dll and .h file in the same directory as grace is in so that the compiler can find them. 

*** How to Run this software ***
Two example scripts were provided, cubic.txt and stdprob4.txt. The first one simulates a cubic magnetic nano particle relaxed under zero field, and the second one is on micromagnetic standard problem #4, field 1. 

Drag and drop <inputFileName>.txt (e.g. cube.txt) to grace.exe and the simulation will begin automatically.

*** What to do while the software is running ***
The program window will show calculation progress (e.g. 10%, 99%). This may take minutes or hours depending on your calculation size.

*** What to do after the program finishes ***
Double click <inputFileName>.plt (e.g. cube.plt) to view visualization of simulated magnetic dynamics. The change of average magnetizations (Mx, My, Mz) is saved in <inputFileName>.mvst.txt. 

Detailed snapshots of the magnetization configuration are saved in <inputFilename>.data.txt. These files may be huge, but you can use professional text editors (Notepad++ etc.) to open them.

*** Questions? ***
Questions or comments may be addressed to zhu@sting.graceland.edu.
