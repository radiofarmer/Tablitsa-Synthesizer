Tablitsa Synthesizer (In Development)
==========================
A polyphonic wavetable synthesizer inspired by chemistry, built with the the IPlug2 framework, a fork of Cockos WDL by Oli Larkin. Available as a stand-alone app or as a VST3 plugin.

## Features
The final version will include
  * Two wavetable oscillators with adjustable table position ("Valency"), phase distortion ("Excitation"), and bass boost and saturation ("Mass")
  * 118 wavetables with properties inspired by chemical elements
  * A filter for each oscillator, using Hal Chamberlain's State-Variable Filter (lowpass, highpass, bandpass, allpass) and a model of the Moog Ladder Filter (lowpass, highpass, bandpass), as well as a comb filter with adjustable feedback and feedforward coefficients and adjustable delay line length
  * Three envelopes, two LFOs, a 16-step sequencer, and keytrack, velocity, and trigger-random modulation for all continuously-variable parameters
  * Phase and ring modulation for each voice
  * Per-voice waveshaping distortion
  * Monophonic mode with adjustable portamento
  * Master effects including delay, phaser, and sample-and-hold
  
Features as of January 14, 2021:
  * Both oscillators are functional
  * Wavetables from Hydrogen to Gold (in order of atomic number) have been generated. Selecting higher-atomic-number wavetables will have no effect, but will not crash the plugin.
  * All envelopes are fully-functional; LFOs are mostly functional; Sequencer control is not functional.
  * State-Variable and Moog filters are functional; comb filter has not yet been added
  
## Compatibility
*built binaries are currently not available*

The stand-alone app has only been tested on Windows 10, and the VST3 plugin only in Cockos REAPER, but it should work on other platforms. The release version uses the AVX2 instruction set and will only work on new-ish processors. A maximum-compatibility version without vector instructions will eventually be included in the release.

## Cloning and Building
You are free to build upon this project according to the license. (Though fair warning, the code is kind of a mess as this is my first real C++ project.) In order to build the solution, first install IPlug2 and its third-party dependencies per the instructions on the [Wiki](https://github.com/iPlug2/iPlug2/wiki). You will also need to download [Agner Fog's `vectorclass` library](https://github.com/vcoda/vectorclass). Place these libraries in a folder within the *Tablitsa* directory called "dependencies". The include directory is specified as `$(ProjectDir)..\vectorclass` - this is where the header files should be located. The WDL FFT function, included with iPlug2, is also used (see top of `Wavetable.h`). Both `WDL\fft.h` and `WDL\fft.c` should be included in the projects already, but if this somehow doesn't work, you can find them in the `iPlug2\WDL\` directory.
