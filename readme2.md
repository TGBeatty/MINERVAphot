# NIRCam Grism Time Series Simulation

Contibutors: Thomas Beatty [tgbeatty@arizona.edu](mailto:tgbeatty@arizona.edu), Everett Schlawin [eas342@email.arizona.edu](mailto:eas342@email.arizona.edu), Bryan Hilbert, Arsh Nadkarni

## About

This directory containts simulated NIRCam grism time series data. It uses the [MIRAGE tool](https://mirage-data-simulator.readthedocs.io) to simulate time series. This includes detector effects (like 1/f noise, amplifier offsets, RC pixels, etc.) from ground-based testing. MIRAGE adds photon noise and cosmic rays. We (the NIRCam instrument team) have modified MIRAGE to include more realizations of the detector noise than in the base version (35 realizations, vs. 5 in base MIRAGE). Note that MIRAGE does not include flat fielding effects (as of late April, 2021) -- though the pipeline also does not currently do a flat field correction (the flat reference is blank until we get on-orbit). 

Here, a WASP-43 b-like planet with artificial spectral features is simulated from the ERS Simulated Spectra Team. The observing mode is as follows:

* GRISMR+F322W2 pupil and filter
* RAPID [readout mode](https://jwst-docs.stsci.edu/near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-readout-patterns)
* 19 Groups per integrations
* 1287 integrations
* 1 Exposure
* 4 Output amplifiers

## Products

These data were calibrated using a development version of the JWST reduction pipeline that was cloned from the [pipeline GitHub](https://github.com/spacetelescope/jwst) on April 14, 2021. This dev version included critical fixes to the reference pixel correction routines that (as of late April 2021) are not currently in a release version of the pipeline.

The calibrated files have been separated into seperate directories containing the [Stage 1 and Stage 2 outputs](https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/stages.html). The expectation is that observers will generally use the Stage 2 outputs, but currently it may be useful to experiment with the Stage 1 files. [The pipeline documentation provides a detailed description of [what each of the output files types are](https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/product_types.html), as well as the exact steps used in [Stage 1](https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_detector1.html#calwebb-detector1) and [Stage 2](https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_spec2.html#calwebb-spec2) processing for TSO data.

Note that for these data we have modified the default STScI pipeline to include a realistic dark frame from ground-testing (currently the NIRCam dark is blank), and changed the detection threshold in the "jump" step of Stage 1 to be 15 (up from the default of 3). 

Briefly:

1. Uncalibrated
  * `Uncalibrated/*uncal.fits` - raw uncalibrated images in measured counts (DN). This contains the counts in every group of every integration up the ramp. Each segment is composed of many integrations.

2. Stage 1

 - `Stage1/*rateints.fits` ; rate images (counts/DN per second) output from the JWST pipeline stage 1 processing. This containts the slope fits to the samples in each integration.
 - `intermediate_products` ; these are intermediate products from the JWST pipeline, which can be viewed for diagnostic purposes but are not needed for producing time series.
 - `intermediate_products/rate` ; these contain the average rate over all the integrations in a given segment. This averages over all time information for that segment.

Note that for practicing on smaller file sizes, one may use segment 21. More information about the JWST pipeline may be found at:
[https://jwst-pipeline.readthedocs.io/en/latest/jwst/introduction.html](https://jwst-pipeline.readthedocs.io/en/latest/jwst/introduction.html)

## Model Parameters (Contains Spoilers!)

The parameters from the modeling team are below:

* Star mass: 0.717 Msun
* Star radius: 0.667 Rsun
* Star temp: 4520 K
* Star is a blackbody
* Planet mass: 0.2 Mjup
* Planet radius: 0.946 Rjup at 100 bars reference pressure (equivalent to 1.036 Rjup at 0.1 bars reference pressure, assuming an isothermal atmosphere at 1400 K)
* CO mixing ratio of 3.1e-4.
* Transmission case is an isothermal atmosphere at 1400 K.
* Emission cases use the TP profiles provided by Patricio.
* Rjup: 6.6854e7 m (Prsa et al., 2016)
* Mjup: 1.89812e27 kg (Prsa et al., 2016)
*
