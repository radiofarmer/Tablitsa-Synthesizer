Tablitsa Synthesizer (In Development)
==========================
A polyphonic wavetable synthesizer inspired by chemistry, built with the the IPlug2 framework, a fork of Cockos WDL by Oli Larkin. Available as a stand-alone app or as a VST3 plugin.

## Features
The final version will include
  * Two wavetable oscillators with adjustable table position/timbre, phase distortion, and formant shift. The wavetables are band-limited by 8x oversampled, octave-spaced mipmaps and by passing the 2x oversampled oscillator output through a 12th-order IIR filter.
  * 118 wavetables with properties inspired by chemical elements: Each group (column) on the table has a characteristic spectral shape that changes slightly with period (row), while different sets of harmonic and inharmonic frequencies are added to positions in the wavetable representing the most common oxidation states for a given element. (e.g. Carbon's wavetable has 9 different timbres, corresponding to oxidation states -4 to +4.) 
  * Two filters with adjustable routing from the oscillators, using Hal Chamberlin's State-Variable Filter (lowpass, highpass, bandpass, allpass) and a model of the Moog Ladder Filter (lowpass, highpass, bandpass, with -12dB/Oct. and -24dB/Oct. slopes), as well as a comb filter with adjustable feedback and feedforward coefficients and adjustable delay line length
  * Three envelopes, two LFOs, a 16-step sequencer, and keytrack, velocity, and trigger-random modulation for all continuously-variable parameters
  * Phase and ring modulation for each voice
  * Monophonic mode with portamento
  * Up to 8 unison voices, which can be detuned up to an octave and distributed across several common chords
  * Modulatable voice effects including sample-and-hold, waveshaping and distortion, and compression
  * Master effects including delay, EQ, and reverb
  
## Compatibility
*See Releases in sidebar for binaries and installation packages*

The stand-alone app is built for Windows has only been tested on Windows 10---iPlug2 apps can be compiled for MacOS but must be compiled in that environment, and I currently don't have access to a machine running MacOS. The VST3 build has only been tested in REAPER but should work in any other DAW which supports the VST3 standard. All release versions use the AVX2 instruction set and will only work on new-ish processors. A maximum-compatibility version without vector instructions will eventually be included in the release. Wavetables are expected to be located in the folder `%AppData%\Tablitsa\wavetables\` (i.e. in the the `AppData\Roaming` folder) in release builds. (For debug builds, their location within the repository is used.) If this doesn't work on your machine, see `Wavetable.h` for the IO functions.

## Cloning and Building
You are free to build upon this project according to the MIT license. (Though fair warning, the code is kind of a mess as this is my first real C++ project.) In order to build the solution, first install *my fork of IPlug2* and the third-party dependencies for IPlug2 per the instructions on the [Wiki](https://github.com/iPlug2/iPlug2/wiki). You will also need to download [Agner Fog's `vectorclass` library](https://github.com/vcoda/vectorclass). Place these libraries in a folder within the *Tablitsa* directory called "dependencies". The include directory is specified as `$(ProjectDir)..\dependencies\vectorclass` - this is where the header files should be located. The WDL FFT function, included with iPlug2, is also used (see top of `Wavetable.h`). Both `WDL\fft.h` and `WDL\fft.c` should be included in the projects already, but if this somehow doesn't work, you can find them in the `iPlug2\WDL\` directory. A library of my own implementations of common DSP systems is also included in `$(ProjectDir)..\dependencies\radiofarmer`.

For tips on modifying the application to better suit your needs (e.g. making new wavetables or adjusting quality-performance tradeoffs) see the wiki (still in progress). The Python script used to generate the wavetables, including an explanation of the `.wt` format and the methodology used to generate and read wavetables, is available in my WavetableEditor repository.

## Support Tablitsa!
If you enjoy using Tablitsa or its source code and would like to support me, you can donate via PayPal:

[![paypal](https://www.paypalobjects.com/en_US/i/btn/btn_donateCC_LG.gif)](https://www.paypal.com/donate?business=NUFQNKX9ET4VL&currency_code=USD)
