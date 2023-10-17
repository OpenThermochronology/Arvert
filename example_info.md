## Arvert 7 Examples

I've provided four different working examples for Arvert 7

These are all sample suites I've worked with before, but I have NOT made an attempt to demonstrate a complete analysis using the model - these are just intended to help you get Arvert running.

For all four examples, the crs.in file has been reset to begin a new model. For additional runs, you can look at the modelinfo file to see what was changed for that run.

For each exampl, I ran an initial model with about 2500 CRS iterations, and then a number of post-processing iterations for the second run. You can open the modelinfo file for each run to see the record of all the parameters that were used.

### Examples
---

* **One sample with one <sup>40</sup>Ar/<sup>39</sup>Ar age spectrum**

	_Kohistan_one_MDD_

	This is just a classic single MDD sample from the NW Himalaya in Pakistan. The sample has a mild domain structure and converges well. There are no mineral data. 2500 CRS iterations for the first run, and for the second run, restarted from the first, there are another 2500 CRS iterations followed by 5000 post-processing iterations.

* **Multiple samples blending an age spectrum and mineral ages along a vertical profile**

	*Kohistan_MDD_minerals*

	This includes the MDD data from the simple example (above) but adds some zircon and apatite U-He data for eight samples spread across almost 3000 meters of steep relief. Thus, these were modeled with T-offsets (in these examples I let these be part of the model, because the extreme relief might have led isotherms to be warped, reducing the T-offsets compared to assuming a normal thermal gradient). 2500 CRS iterations for the first run, and for the second run, restarted from the first, there are another 2500 CRS iterations followed by 5000 post-processing iterations.

* **Multiple samples from a borehole each having a <sup>40</sup>Ar/<sup>39</sup>Ar age spectrum**

	*URL_only_multiple_MDD*

	These samples are from a ~1 km borehole into crystalline rocks of the Canadian shield. There are six K-feldspar MDD age spectra and no mineral data. A mild monotonic T-offset would be expected for these samples, given that they are from a borehole. 2500 CRS iterations for the first run, and for the second run, restarted from the first, there are no additional CRS iterations followed by 10000 post-processing iterations.

* **Multiple samples having no age spectrum but several mineral ages each, all from surface samples**

	*sturrock_only_minerals*

	This is very large set of apatite U-He data from the Sturrock et al. paper (2021; https://doi.org/10.1029/2020GC009567), spread across a number of surface samples. Thus, all T-offsets were fixed to be 0.0. 2500 CRS iterations for the first run, and for the second run, restarted from the first, there are no additional CRS iterations followed by 10000 post-processing iterations.
