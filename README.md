[![DOI](https://zenodo.org/badge/55851714.svg)](https://zenodo.org/badge/latestdoi/55851714)

# RIR-Generator

The image method, proposed by Allen and Berkley in 1979 [1], is probably one of the most frequently used methods in the acoustic signal processing community to create synthetic room impulse responses. 
A mex-function, which can be used in MATLAB, was developed to generate multi-channel room impulse responses using the image method. 
This function enables the user to control the reflection order, room dimension, and microphone directivity. 

This repository includes a tutorial, MATLAB examples, and the source code of the mex-function.

More information can be found [here](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator).

# Build the MEX function

Below is a concise, platform-specific guide. First, follow the setup for your OS, then the common compile command.

## Platform setup

### Windows (x64)
1. Install a supported MATLAB version (R2018b or later recommended).
2. Install a supported C++ compiler, e.g.:
	- Microsoft Visual Studio 2019/2022 with "Desktop development with C++", or
	- Microsoft C++ Build Tools (standalone)
3. In MATLAB, select the compiler:
	```matlab
	mex -setup C++
	```
4. Compile using the common compile command below.

### macOS Apple Silicon (ARM64)
1. Install a native Apple Silicon MATLAB (R2023b or later).
2. Install Xcode Command Line Tools in Terminal:
	```zsh
	xcode-select --install
	```
3. In MATLAB, configure the compiler:
	```matlab
	mex -setup C++
	```
4. Compile using the common compile command below.

### macOS Intel (x86_64) or MATLAB under Rosetta
1. Use an Intel MATLAB (e.g., R2020b or later) or run MATLAB under Rosetta on Apple Silicon.
2. Install Xcode Command Line Tools:
	```zsh
	xcode-select --install
	```
3. In MATLAB, configure the compiler:
	```matlab
	mex -setup C++
	```
4. Compile using the common compile command below.

### Linux (x86_64)
1. Install GCC toolchain (e.g., build-essential on Debian/Ubuntu).
2. In MATLAB, select GCC:
	```matlab
	mex -setup C++
	```
3. Compile using the common compile command below.

## Common compile command (all platforms)

Run this in MATLAB from the repository folder:

```matlab
mex rir_generator.cpp rir_generator_core.cpp
```

Resulting MEX extension by platform:
- Windows (x64): `.mexw64`
- macOS Apple Silicon (ARM64): `.mexmaca64`
- macOS Intel or MATLAB under Rosetta (x86_64): `.mexmaci64`
- Linux (x86_64): `.mexa64`

# References

1. J.B. Allen and D.A. Berkley, "Image method for efficiently simulating small-room acoustics," Journal Acoustic Society of America, 65(4), April 1979, p 943.
