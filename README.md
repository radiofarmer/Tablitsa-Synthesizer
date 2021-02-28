Tablitsa Synthesizer (In Development)
==========================
A polyphonic wavetable synthesizer inspired by chemistry, built with the the IPlug2 framework, a fork of Cockos WDL by Oli Larkin. Available as a stand-alone app or as a VST3 plugin.

## Features
The final version will include
  * Two wavetable oscillators with adjustable table position/timbre, phase distortion, and formant shift
  * 118 wavetables with properties inspired by chemical elements: Each group (column) on the table has a characteristic spectral shape that changes slightly with period (row), while different sets of harmonic and inharmonic frequencies are added to positions in the wavetable representing the most common oxidation states for a given element. (e.g. Carbon's wavetable has 9 different timbres, corresponding to oxidation states -4 to +4.) 
  * A filter for each oscillator, using Hal Chamberlain's State-Variable Filter (lowpass, highpass, bandpass, allpass) and a model of the Moog Ladder Filter (lowpass, highpass, bandpass, with -12dB/Oct. and -24dB/Oct. slopes), as well as a comb filter with adjustable feedback and feedforward coefficients and adjustable delay line length
  * Three envelopes, two LFOs, a 16-step sequencer, and keytrack, velocity, and trigger-random modulation for all continuously-variable parameters
  * Phase and ring modulation for each voice
  * Monophonic mode with portamento
  * Up to 8 unison voices, which can be detuned up to an octave and distributed across several common chords
  * Master effects including delay and sample-and-hold
  
## Compatibility
*built binaries are currently not available in the repository*

The stand-alone app is built for Windows has only been tested on Windows 10---iPlug2 apps can be compiled for MacOS but must be compiled in that environment, and I currently don't have access to a machine running MacOS. The VST3 build has only been tested in REAPER but should work in any other DAW which supports the VST3 standard. All release versions use the AVX2 instruction set and will only work on new-ish processors. A maximum-compatibility version without vector instructions will eventually be included in the release. Wavetables are expected to be located in the folder `%AppData%\Tablitsa\wavetables\` (i.e. in the the `AppData\Roaming` folder) in release builds. (For debug builds, their location within the repository is used.) If this doesn't work on your machine, see `Wavetable.h` for the IO functions.

## Cloning and Building
You are free to build upon this project according to the license. (Though fair warning, the code is kind of a mess as this is my first real C++ project.) In order to build the solution, first install *my fork of IPlug2* and the third-party dependencies for IPlug2 per the instructions on the [Wiki](https://github.com/iPlug2/iPlug2/wiki). You will also need to download [Agner Fog's `vectorclass` library](https://github.com/vcoda/vectorclass). Place these libraries in a folder within the *Tablitsa* directory called "dependencies". The include directory is specified as `$(ProjectDir)..\dependencies\vectorclass` - this is where the header files should be located. The WDL FFT function, included with iPlug2, is also used (see top of `Wavetable.h`). Both `WDL\fft.h` and `WDL\fft.c` should be included in the projects already, but if this somehow doesn't work, you can find them in the `iPlug2\WDL\` directory.

For tips on modifying the application to better suit your needs (e.g. making new wavetables or adjusting quality-performance tradeoffs) see the wiki (still in progress). The Python script used to generate the wavetables, including an explanation of the `.wt` format and the methodology used to generate and read wavetables, is available in my WavetableEditor repository.
